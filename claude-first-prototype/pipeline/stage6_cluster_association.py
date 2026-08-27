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
    "cluster_haplotype_count", "carrier_rate_with_cluster",
    "carrier_rate_without_cluster", "population_adjusted_r", "p_value",
    "q_value", "informative_populations", "directional_consistency",
]


def dosage_table(sv: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    result = sv[METADATA_COLUMNS].copy()
    for sample in samples:
        alternate, called = genotype_counts(sv[sample])
        result[sample] = np.where(called > 0, alternate, np.nan)
    return result


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


def empirical_p_value(
    x: np.ndarray,
    y: np.ndarray,
    populations: np.ndarray,
    observed: float,
    permutations: int,
    rng: np.random.Generator,
) -> float:
    if not np.isfinite(observed):
        return 1.0
    exceedances = 0
    shuffled = y.copy()
    for _ in range(permutations):
        for population in np.unique(populations):
            use = np.flatnonzero(populations == population)
            shuffled[use] = rng.permutation(y[use])
        permuted = correlation(x, shuffled, populations)
        exceedances += np.isfinite(permuted) and abs(permuted) >= abs(observed)
    return (exceedances + 1) / (permutations + 1)


def association_table(
    sv: pd.DataFrame,
    sv_blocks: pd.DataFrame,
    memberships: pd.DataFrame,
    metadata: pd.DataFrame,
    permutations: int,
    seed: int,
    min_cluster_haplotypes: int,
    min_population_samples: int,
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
            for cluster_id in cluster_dosage.columns:
                cluster_samples = cluster_dosage.index[cluster_dosage[cluster_id] > 0]
                if int(cluster_dosage[cluster_id].sum()) < min_cluster_haplotypes:
                    continue
                called_samples = [
                    sample for sample in cluster_dosage.index
                    if sample in samples and pd.notna(variant[sample])
                ]
                if len(called_samples) < 4:
                    continue
                x = cluster_dosage.loc[called_samples, cluster_id].to_numpy(dtype=float)
                y = variant[called_samples].to_numpy(dtype=float)
                populations = np.array([population_by_sample[sample] for sample in called_samples])
                if len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
                    continue
                adjusted_r = correlation(x, y, populations)
                p_value = empirical_p_value(
                    x, y, populations, adjusted_r, permutations, rng
                )
                informative, consistency = population_consistency(
                    x, y, populations, adjusted_r, min_population_samples
                )
                with_cluster = y[x > 0] > 0
                without_cluster = y[x == 0] > 0
                rows.append(
                    {
                        **{column: variant[column] for column in METADATA_COLUMNS},
                        "haploblock_id": block_id,
                        "cluster_id": cluster_id,
                        "n_called": len(called_samples),
                        "n_cluster_carriers": len(cluster_samples.intersection(called_samples)),
                        "cluster_haplotype_count": int(x.sum()),
                        "carrier_rate_with_cluster": float(with_cluster.mean()),
                        "carrier_rate_without_cluster": float(without_cluster.mean()),
                        "population_adjusted_r": adjusted_r,
                        "p_value": p_value,
                        "informative_populations": informative,
                        "directional_consistency": consistency,
                    }
                )

    if not rows:
        return pd.DataFrame(columns=ASSOCIATION_COLUMNS)
    result = pd.DataFrame(rows)
    result["q_value"] = benjamini_hochberg(result["p_value"].to_numpy())
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
            "population_adjusted_r", "q_value", "informative_populations",
            "directional_consistency",
        ])
    best = (
        associations.assign(abs_r=associations["population_adjusted_r"].abs())
        .sort_values(["q_value", "abs_r"], ascending=[True, False])
        .drop_duplicates(["sv_record_id", "haploblock_id"])
        .copy()
    )
    significant = (best["q_value"] < q_threshold) & (best["abs_r"] >= min_abs_r)
    portable = (
        significant
        & (best["informative_populations"] >= 2)
        & (best["directional_consistency"] >= 0.75)
    )
    dependent = significant & (best["informative_populations"] >= 2) & ~portable
    best["association_pattern"] = "no_detected_cluster_signal"
    best.loc[significant, "association_pattern"] = "cluster_associated"
    best.loc[dependent, "association_pattern"] = "population_dependent_association"
    best.loc[portable, "association_pattern"] = "portable_cluster_tag"
    return best.rename(columns={"cluster_id": "best_cluster_id"}).drop(columns="abs_r")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("stage6_output"))
    parser.add_argument("--permutations", type=int, default=1000)
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
            )
        )
    associations = pd.concat(all_associations, ignore_index=True)
    if not associations.empty:
        associations["q_value"] = benjamini_hochberg(associations["p_value"].to_numpy())
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
