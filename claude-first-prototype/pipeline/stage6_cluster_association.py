"""Stage 6: test local SV-cluster association beyond population membership.

For each overlapping SV-haploblock pair, local cluster dosage is compared with
SV dosage after both have been centered within population. Empirical p-values
shuffle SV dosages within population, preserving the population allele-frequency
pattern that made a simple population-cluster agreement test difficult to
interpret biologically.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from stage4_classify_af import genotype_counts, resolve_path
from stage5_type_enrichment import benjamini_hochberg
from sv_contract import METADATA_COLUMNS


ASSOCIATION_COLUMNS = [
    *METADATA_COLUMNS,
    "haploblock_id", "cluster_id", "n_called", "n_cluster_carriers",
    "n_samples_with_cluster", "n_sv_carriers_with_cluster",
    "n_sv_noncarriers_with_cluster", "cluster_haplotype_count",
    "carrier_rate_with_cluster", "carrier_rate_without_cluster",
    "carrier_rate_difference", "association_direction",
    "population_adjusted_r", "p_value", "q_value", "permutations_used",
    "informative_populations", "directional_consistency",
]


def dosage_table(sv: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    result = sv[METADATA_COLUMNS].copy()
    dosage_columns = {}
    for sample in samples:
        alternate, called = genotype_counts(sv[sample])
        dosage_columns[sample] = np.where(called > 0, alternate, np.nan)
    return pd.concat([result, pd.DataFrame(dosage_columns)], axis=1)


def centered(values: np.ndarray, populations: np.ndarray) -> np.ndarray:
    result = values.astype(float).copy()
    for population in np.unique(populations):
        use = populations == population
        result[use] -= result[use].mean()
    return result


def correlation(x: np.ndarray, y: np.ndarray, populations: np.ndarray) -> float:
    x_residual = centered(x, populations)
    y_residual = centered(y, populations)
    denominator = np.sqrt(np.sum(x_residual**2) * np.sum(y_residual**2))
    return float(np.sum(x_residual * y_residual) / denominator) if denominator else np.nan


def population_consistency(
    x: np.ndarray,
    y: np.ndarray,
    populations: np.ndarray,
    overall_r: float,
    min_population_samples: int,
) -> tuple[int, float]:
    effects = []
    for population in np.unique(populations):
        use = populations == population
        if use.sum() < min_population_samples or len(np.unique(x[use])) < 2:
            continue
        with_cluster = y[use & (x > 0)]
        without_cluster = y[use & (x == 0)]
        if len(with_cluster) and len(without_cluster):
            effects.append(with_cluster.mean() - without_cluster.mean())
    if not effects or not np.isfinite(overall_r) or overall_r == 0:
        return len(effects), np.nan
    consistency = np.mean(np.sign(effects) == np.sign(overall_r))
    return len(effects), float(consistency)


def max_statistic_p_values(
    cluster_vectors: dict[str, np.ndarray],
    y: np.ndarray,
    populations: np.ndarray,
    observed: dict[str, float],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[dict[str, float], int]:
    """Family-wise p-values using the largest cluster statistic per permutation."""
    finite_observed = [abs(value) for value in observed.values() if np.isfinite(value)]
    if not finite_observed:
        return {cluster_id: 1.0 for cluster_id in observed}, 0
    maxima = np.zeros(permutations, dtype=float)
    cluster_matrix = np.column_stack(list(cluster_vectors.values())).astype(float)
    y_residual = centered(y, populations)
    for population in np.unique(populations):
        use = populations == population
        cluster_matrix[use] -= cluster_matrix[use].mean(axis=0)
    cluster_norms = np.sqrt(np.sum(cluster_matrix**2, axis=0))
    population_indices = [
        np.flatnonzero(populations == population) for population in np.unique(populations)
    ]
    chunk_size = 1000
    for start in range(0, permutations, chunk_size):
        width = min(chunk_size, permutations - start)
        shuffled = np.empty((len(y), width), dtype=float)
        for indices in population_indices:
            values = np.repeat(y_residual[indices, None], width, axis=1)
            shuffled[indices] = rng.permuted(values, axis=0)
        y_norms = np.sqrt(np.sum(shuffled**2, axis=0))
        denominators = cluster_norms[:, None] * y_norms[None, :]
        correlations = np.divide(
            cluster_matrix.T @ shuffled,
            denominators,
            out=np.full_like(denominators, np.nan),
            where=denominators > 0,
        )
        maxima[start:start + width] = np.nanmax(np.abs(correlations), axis=0)
    return {
        cluster_id: (
            float((np.sum(maxima >= abs(value)) + 1) / (permutations + 1))
            if np.isfinite(value) else 1.0
        )
        for cluster_id, value in observed.items()
    }, permutations


def assign_pair_q_values(associations: pd.DataFrame) -> pd.DataFrame:
    """Correct one maximum-statistic p-value per SV–block pair."""
    if associations.empty:
        return associations
    keys = ["sv_record_id", "haploblock_id"]
    pair_tests = associations.groupby(keys, as_index=False)["p_value"].min()
    pair_tests["q_value"] = benjamini_hochberg(pair_tests["p_value"].to_numpy())
    return associations.drop(columns="q_value", errors="ignore").merge(
        pair_tests[keys + ["q_value"]], on=keys, how="left"
    )


def association_table(
    sv: pd.DataFrame,
    sv_blocks: pd.DataFrame,
    memberships: pd.DataFrame,
    metadata: pd.DataFrame,
    permutations: int,
    seed: int,
    min_cluster_haplotypes: int,
    min_population_samples: int,
    refinement_permutations: int | None = None,
    refinement_p_threshold: float = 0.01,
) -> pd.DataFrame:
    samples = [
        sample for sample in metadata["sample_id"].astype(str)
        if sample in sv.columns[len(METADATA_COLUMNS):]
    ]
    population_by_sample = metadata.set_index("sample_id")["population"].astype(str).to_dict()
    variant_key = ["sv_record_id"]
    dosages = dosage_table(sv, samples).set_index(variant_key, drop=False)
    memberships = memberships[memberships["sample_id"].astype(str).isin(samples)].copy()
    memberships["sample_id"] = memberships["sample_id"].astype(str)
    rng = np.random.default_rng(seed)
    rows = []

    for block_id, block_variants in sv_blocks.drop_duplicates(
        [*variant_key, "haploblock_id"]
    ).groupby("haploblock_id", sort=False):
        block_memberships = memberships[memberships["haploblock_id"] == block_id]
        if block_memberships.empty:
            continue
        cluster_dosage = (
            block_memberships.assign(cluster_haplotypes=1)
            .pivot_table(
                index="sample_id", columns="cluster_id", values="cluster_haplotypes",
                aggfunc="sum", fill_value=0,
            )
        )
        for _, block_variant in block_variants.iterrows():
            key = block_variant["sv_record_id"]
            if key not in dosages.index:
                continue
            variant = dosages.loc[key]
            called_samples = [
                sample for sample in cluster_dosage.index
                if sample in samples and pd.notna(variant[sample])
            ]
            if len(called_samples) < 4:
                continue
            y = variant[called_samples].to_numpy(dtype=float)
            populations = np.array([population_by_sample[sample] for sample in called_samples])
            if len(np.unique(y)) < 2:
                continue
            cluster_vectors = {}
            for cluster_id in cluster_dosage.columns:
                if int(cluster_dosage[cluster_id].sum()) < min_cluster_haplotypes:
                    continue
                x = cluster_dosage.loc[called_samples, cluster_id].to_numpy(dtype=float)
                if len(np.unique(x)) < 2:
                    continue
                cluster_vectors[cluster_id] = x
            observed = {
                cluster_id: correlation(x, y, populations)
                for cluster_id, x in cluster_vectors.items()
            }
            if not observed:
                continue
            p_values, permutations_used = max_statistic_p_values(
                cluster_vectors, y, populations, observed, permutations, rng
            )
            if (
                refinement_permutations is not None
                and refinement_permutations > permutations
                and min(p_values.values()) <= refinement_p_threshold
            ):
                # The refinement draws are independent and replace, rather than
                # reuse, the screening p-value selected in the first pass.
                p_values, permutations_used = max_statistic_p_values(
                    cluster_vectors, y, populations, observed,
                    refinement_permutations, rng,
                )
            for cluster_id, x in cluster_vectors.items():
                adjusted_r = observed[cluster_id]
                informative, consistency = population_consistency(
                    x, y, populations, adjusted_r, min_population_samples
                )
                with_cluster = y[x > 0] > 0
                without_cluster = y[x == 0] > 0
                rate_with = float(with_cluster.mean())
                rate_without = float(without_cluster.mean())
                rate_difference = rate_with - rate_without
                if rate_difference > 0:
                    direction = "carrier_enriched"
                elif rate_difference < 0:
                    direction = "carrier_depleted"
                else:
                    direction = "no_difference"
                rows.append(
                    {
                        **{column: variant[column] for column in METADATA_COLUMNS},
                        "haploblock_id": block_id,
                        "cluster_id": cluster_id,
                        "n_called": len(called_samples),
                        # Retained for compatibility; this historically counted called
                        # samples carrying the cluster, not SV carriers.
                        "n_cluster_carriers": int((x > 0).sum()),
                        "n_samples_with_cluster": int((x > 0).sum()),
                        "n_sv_carriers_with_cluster": int(with_cluster.sum()),
                        "n_sv_noncarriers_with_cluster": int((~with_cluster).sum()),
                        "cluster_haplotype_count": int(x.sum()),
                        "carrier_rate_with_cluster": rate_with,
                        "carrier_rate_without_cluster": rate_without,
                        "carrier_rate_difference": rate_difference,
                        "association_direction": direction,
                        "population_adjusted_r": adjusted_r,
                        "p_value": p_values[cluster_id],
                        "permutations_used": permutations_used,
                        "informative_populations": informative,
                        "directional_consistency": consistency,
                    }
                )

    if not rows:
        return pd.DataFrame(columns=ASSOCIATION_COLUMNS)
    result = assign_pair_q_values(pd.DataFrame(rows))
    return result[ASSOCIATION_COLUMNS].sort_values(
        ["q_value", "population_adjusted_r"], ascending=[True, False]
    ).reset_index(drop=True)


def summarize_associations(
    associations: pd.DataFrame,
    q_threshold: float,
    min_abs_r: float,
) -> pd.DataFrame:
    if associations.empty:
        return pd.DataFrame(columns=[
            "sv_id", "haploblock_id", "best_cluster_id", "association_pattern",
            "best_enriched_cluster_id", "best_depleted_cluster_id",
            "association_direction", "carrier_rate_difference",
            "population_adjusted_r", "q_value", "informative_populations",
            "directional_consistency",
        ])
    keys = ["sv_record_id", "haploblock_id"]
    ranked = associations.assign(abs_r=associations["population_adjusted_r"].abs())
    best = (
        ranked.sort_values(["q_value", "p_value", "abs_r"], ascending=[True, True, False])
        .drop_duplicates(["sv_record_id", "haploblock_id"])
        .copy()
    )
    enriched_ids = (
        ranked[ranked["association_direction"] == "carrier_enriched"]
        .sort_values(["p_value", "population_adjusted_r"], ascending=[True, False])
        .drop_duplicates(keys)[keys + ["cluster_id"]]
        .rename(columns={"cluster_id": "best_enriched_cluster_id"})
    )
    depleted_ids = (
        ranked[ranked["association_direction"] == "carrier_depleted"]
        .sort_values(["p_value", "population_adjusted_r"], ascending=[True, True])
        .drop_duplicates(keys)[keys + ["cluster_id"]]
        .rename(columns={"cluster_id": "best_depleted_cluster_id"})
    )
    best = best.merge(enriched_ids, on=keys, how="left").merge(
        depleted_ids, on=keys, how="left"
    )
    significant = (best["q_value"] < q_threshold) & (best["abs_r"] >= min_abs_r)
    enriched = significant & best["association_direction"].eq("carrier_enriched")
    portable = (
        enriched
        & (best["informative_populations"] >= 2)
        & (best["directional_consistency"] >= 0.75)
    )
    dependent = enriched & (best["informative_populations"] >= 2) & ~portable
    best["association_pattern"] = "no_detected_cluster_signal"
    best.loc[significant & ~enriched, "association_pattern"] = "cluster_exclusion_signal"
    best.loc[enriched, "association_pattern"] = "cluster_associated"
    best.loc[dependent, "association_pattern"] = "population_dependent_association"
    best.loc[portable, "association_pattern"] = "cross_population_consistent_tag_candidate"
    return best.rename(columns={"cluster_id": "best_cluster_id"}).drop(
        columns="abs_r"
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("stage6_output"))
    parser.add_argument(
        "--permutations", type=int, default=10000,
        help="Screening maximum-statistic permutations per SV–block.",
    )
    parser.add_argument(
        "--refinement-permutations", type=int, default=1000000,
        help="Independent permutations for pairs passing the refinement threshold.",
    )
    parser.add_argument("--refinement-p-threshold", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--min-cluster-haplotypes", type=int, default=6)
    parser.add_argument("--min-population-samples", type=int, default=4)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--min-abs-r", type=float, default=0.3)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent
    metadata = pd.read_csv(resolve_path(config["paths"]["sample_metadata"], config_dir), sep="\t")
    all_associations = []
    for chrom, genotype_path in config["paths"]["sv_genotypes"].items():
        sv = pd.read_csv(resolve_path(genotype_path, config_dir), sep="\t")
        sv_blocks = pd.read_csv(
            resolve_path(config["paths"]["sv_block_summary"][chrom], config_dir), sep="\t"
        )
        memberships = pd.read_csv(
            resolve_path(config["paths"]["cluster_memberships"][chrom], config_dir), sep="\t"
        )
        all_associations.append(
            association_table(
                sv, sv_blocks, memberships, metadata, args.permutations, args.seed,
                args.min_cluster_haplotypes, args.min_population_samples,
                args.refinement_permutations, args.refinement_p_threshold,
            )
        )
    associations = pd.concat(all_associations, ignore_index=True)
    if not associations.empty:
        associations = assign_pair_q_values(associations)
    summary = summarize_associations(associations, args.q_threshold, args.min_abs_r)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    association_path = args.out_dir / "sv_cluster_associations.tsv"
    summary_path = args.out_dir / "sv_cluster_summary.tsv"
    associations.to_csv(association_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    output_config = dict(config)
    output_config["paths"] = dict(config["paths"])
    output_config["paths"].update({
        "sv_cluster_associations": str(association_path.resolve()),
        "sv_cluster_summary": str(summary_path.resolve()),
    })
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(output_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
