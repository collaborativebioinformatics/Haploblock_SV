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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from stage4_classify_af import resolve_path
from stage5_type_enrichment import benjamini_hochberg
from sv_contract import METADATA_COLUMNS
from sv_dosage import dosage_matrix_from_genotypes


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


def population_adjusted_correlations(
    y: np.ndarray,
    called: np.ndarray,
    cluster_matrix: np.ndarray,
    populations: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Correlate every variant and cluster after centering within population."""
    y_filled = np.where(called, y, 0.0)
    covariance = np.zeros((len(y), cluster_matrix.shape[1]), dtype=float)
    cluster_variance = np.zeros_like(covariance)
    variant_variance = np.zeros(len(y), dtype=float)
    for population in np.unique(populations):
        use = populations == population
        population_called = called[:, use].astype(float)
        population_y = y_filled[:, use]
        population_clusters = cluster_matrix[use]
        counts = population_called.sum(axis=1)
        inverse_counts = np.divide(
            1.0, counts, out=np.zeros_like(counts), where=counts > 0
        )
        y_sum = population_y.sum(axis=1)
        cluster_sum = population_called @ population_clusters
        covariance += (
            population_y @ population_clusters
            - y_sum[:, None] * cluster_sum * inverse_counts[:, None]
        )
        variant_variance += (
            np.square(population_y).sum(axis=1)
            - np.square(y_sum) * inverse_counts
        )
        cluster_variance += (
            population_called @ np.square(population_clusters)
            - np.square(cluster_sum) * inverse_counts[:, None]
        )
    denominator = np.sqrt(cluster_variance * variant_variance[:, None])
    correlations = np.divide(
        covariance,
        denominator,
        out=np.full_like(covariance, np.nan),
        where=denominator > 0,
    )
    return correlations, cluster_variance


def batched_max_statistic_p_values(
    y: np.ndarray,
    called: np.ndarray,
    cluster_matrix: np.ndarray,
    populations: np.ndarray,
    observed: np.ndarray,
    eligible: np.ndarray,
    permutations: int,
    seed: int,
    phase: int,
) -> np.ndarray:
    """Maximum-statistic p-values batched across variants with equal call masks."""
    p_values = np.full_like(observed, np.nan, dtype=float)
    mask_groups: dict[bytes, list[int]] = {}
    for variant_index in np.flatnonzero(eligible.any(axis=1)):
        mask_groups.setdefault(called[variant_index].tobytes(), []).append(variant_index)

    for group_index, variant_indices_list in enumerate(mask_groups.values()):
        variant_indices = np.array(variant_indices_list, dtype=int)
        sample_mask = called[variant_indices[0]]
        group_eligible = eligible[variant_indices[0]]
        cluster_indices = np.flatnonzero(group_eligible)
        group_y = y[np.ix_(variant_indices, sample_mask)].astype(float)
        group_populations = populations[sample_mask]
        group_clusters = cluster_matrix[np.ix_(sample_mask, cluster_indices)].astype(float)
        for population in np.unique(group_populations):
            use = group_populations == population
            group_y[:, use] -= group_y[:, use].mean(axis=1, keepdims=True)
            group_clusters[use] -= group_clusters[use].mean(axis=0, keepdims=True)

        y_norms = np.sqrt(np.square(group_y).sum(axis=1))
        cluster_norms = np.sqrt(np.square(group_clusters).sum(axis=0))
        denominator = y_norms[:, None] * cluster_norms[None, :]
        absolute_observed = np.abs(observed[np.ix_(variant_indices, cluster_indices)])
        exceedances = np.zeros_like(absolute_observed, dtype=np.int64)
        rng = np.random.default_rng(np.random.SeedSequence([seed, phase, group_index]))
        population_indices = [
            np.flatnonzero(group_populations == population)
            for population in np.unique(group_populations)
        ]
        max_elements = 4_000_000
        chunk_size = min(
            1000,
            max(1, max_elements // (len(variant_indices) * len(group_populations))),
        )
        for start in range(0, permutations, chunk_size):
            width = min(chunk_size, permutations - start)
            shuffled = np.empty(
                (len(variant_indices), len(group_populations), width), dtype=float
            )
            for indices in population_indices:
                values = np.broadcast_to(
                    group_y[:, indices, None],
                    (len(variant_indices), len(indices), width),
                )
                shuffled[:, indices] = rng.permuted(values, axis=1)
            correlations = np.divide(
                np.einsum("nc,vnp->vcp", group_clusters, shuffled, optimize=True),
                denominator[:, :, None],
                out=np.full(
                    (len(variant_indices), len(cluster_indices), width), np.nan
                ),
                where=denominator[:, :, None] > 0,
            )
            maxima = np.max(
                np.where(np.isfinite(correlations), np.abs(correlations), -np.inf),
                axis=1,
            )
            exceedances += np.sum(
                maxima[:, None, :] >= absolute_observed[:, :, None], axis=2
            )
        group_p_values = (exceedances + 1) / (permutations + 1)
        group_p_values[~np.isfinite(absolute_observed)] = 1.0
        p_values[np.ix_(variant_indices, cluster_indices)] = group_p_values
    return p_values


def staged_permutation_p_values(
    y: np.ndarray,
    called: np.ndarray,
    cluster_matrix: np.ndarray,
    populations: np.ndarray,
    observed: np.ndarray,
    eligible: np.ndarray,
    permutations: int,
    refinement_permutations: int | None,
    refinement_p_threshold: float,
    seed: int,
    initial_permutations: int | None = None,
    initial_p_threshold: float = 0.1,
    min_refinement_abs_r: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run optional triage, screening, and independent refinement permutations."""
    effect_eligible = eligible.any(axis=1)
    if min_refinement_abs_r is not None:
        maximum_effect = np.max(
            np.where(eligible & np.isfinite(observed), np.abs(observed), -np.inf),
            axis=1,
        )
        effect_eligible &= maximum_effect >= min_refinement_abs_r

    use_initial = (
        initial_permutations is not None
        and 0 < initial_permutations < permutations
    )
    first_permutations = initial_permutations if use_initial else permutations
    p_values = batched_max_statistic_p_values(
        y, called, cluster_matrix, populations, observed, eligible,
        first_permutations, seed, phase=0,
    )
    permutations_used = np.where(eligible.any(axis=1), first_permutations, 0)

    if use_initial:
        row_minimum = np.min(np.where(eligible, p_values, np.inf), axis=1)
        advance = effect_eligible & (row_minimum <= initial_p_threshold)
        if advance.any():
            screen_eligible = eligible & advance[:, None]
            screened = batched_max_statistic_p_values(
                y, called, cluster_matrix, populations, observed, screen_eligible,
                permutations, seed, phase=1,
            )
            p_values[screen_eligible] = screened[screen_eligible]
            permutations_used[advance] = permutations

    if refinement_permutations is not None and refinement_permutations > permutations:
        row_minimum = np.min(np.where(eligible, p_values, np.inf), axis=1)
        refine = effect_eligible & (row_minimum <= refinement_p_threshold)
        if refine.any():
            refinement_eligible = eligible & refine[:, None]
            refined = batched_max_statistic_p_values(
                y, called, cluster_matrix, populations, observed, refinement_eligible,
                refinement_permutations, seed, phase=2,
            )
            p_values[refinement_eligible] = refined[refinement_eligible]
            permutations_used[refine] = refinement_permutations
    return p_values, permutations_used


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


def association_rows_for_block(
    block_id: str,
    block_variants: pd.DataFrame,
    block_memberships: pd.DataFrame,
    sv: pd.DataFrame,
    dosage_matrix: np.ndarray,
    record_positions: dict[str, int],
    sample_positions: dict[str, int],
    population_by_sample: dict[str, str],
    permutations: int,
    seed: int,
    min_cluster_haplotypes: int,
    min_population_samples: int,
    refinement_permutations: int | None,
    refinement_p_threshold: float,
    initial_permutations: int | None,
    initial_p_threshold: float,
    min_refinement_abs_r: float | None,
) -> pd.DataFrame:
    """Test all SVs in one haploblock with an independent random stream."""
    if block_memberships.empty:
        return pd.DataFrame()
    cluster_dosage = (
        block_memberships.assign(cluster_haplotypes=1)
        .pivot_table(
            index="sample_id", columns="cluster_id", values="cluster_haplotypes",
            aggfunc="sum", fill_value=0,
        )
    )
    block_samples = cluster_dosage.index.to_list()
    block_sample_positions = np.array(
        [sample_positions[sample] for sample in block_samples], dtype=int
    )
    block_populations = np.array(
        [population_by_sample[sample] for sample in block_samples]
    )
    cluster_ids = cluster_dosage.columns.to_numpy()
    block_cluster_matrix = cluster_dosage.to_numpy(dtype=float)
    eligible_clusters = block_cluster_matrix.sum(axis=0) >= min_cluster_haplotypes
    cluster_ids = cluster_ids[eligible_clusters]
    block_cluster_matrix = block_cluster_matrix[:, eligible_clusters]
    if not len(cluster_ids):
        return pd.DataFrame()

    variant_positions = np.array([
        record_positions[record_id]
        for record_id in block_variants["sv_record_id"]
        if record_id in record_positions
    ], dtype=int)
    if not len(variant_positions):
        return pd.DataFrame()
    y = dosage_matrix[np.ix_(variant_positions, block_sample_positions)]
    called = np.isfinite(y)
    called_float = called.astype(float)
    y_filled = np.where(called, y, 0.0)
    n_called = called.sum(axis=1)
    y_minimum = np.min(np.where(called, y, np.inf), axis=1)
    y_maximum = np.max(np.where(called, y, -np.inf), axis=1)
    variant_eligible = (n_called >= 4) & (y_minimum < y_maximum)

    observed, _ = population_adjusted_correlations(
        y, called, block_cluster_matrix, block_populations
    )
    cluster_sum = called_float @ block_cluster_matrix
    cluster_sum_squares = called_float @ np.square(block_cluster_matrix)
    global_cluster_variance = cluster_sum_squares - np.divide(
        np.square(cluster_sum), n_called[:, None],
        out=np.zeros_like(cluster_sum), where=n_called[:, None] > 0,
    )
    eligible = variant_eligible[:, None] & (global_cluster_variance > 0)
    if not eligible.any():
        return pd.DataFrame()

    p_values, permutations_used = staged_permutation_p_values(
        y, called, block_cluster_matrix, block_populations, observed, eligible,
        permutations, refinement_permutations, refinement_p_threshold, seed,
        initial_permutations, initial_p_threshold, min_refinement_abs_r,
    )

    cluster_present = (block_cluster_matrix > 0).astype(float)
    sv_carriers = (called & (y > 0)).astype(float)
    n_samples_with_cluster = called_float @ cluster_present
    n_sv_carriers_with_cluster = sv_carriers @ cluster_present
    n_sv_carriers = sv_carriers.sum(axis=1)
    n_samples_without_cluster = n_called[:, None] - n_samples_with_cluster
    n_sv_carriers_without_cluster = (
        n_sv_carriers[:, None] - n_sv_carriers_with_cluster
    )
    carrier_rate_with = np.divide(
        n_sv_carriers_with_cluster, n_samples_with_cluster,
        out=np.full_like(n_samples_with_cluster, np.nan),
        where=n_samples_with_cluster > 0,
    )
    carrier_rate_without = np.divide(
        n_sv_carriers_without_cluster, n_samples_without_cluster,
        out=np.full_like(n_samples_without_cluster, np.nan),
        where=n_samples_without_cluster > 0,
    )
    rate_difference = carrier_rate_with - carrier_rate_without
    direction = np.full(rate_difference.shape, "no_difference", dtype=object)
    direction[rate_difference > 0] = "carrier_enriched"
    direction[rate_difference < 0] = "carrier_depleted"
    direction[~np.isfinite(rate_difference)] = "unavailable_comparator"

    informative_populations = np.zeros_like(observed, dtype=int)
    consistent_populations = np.zeros_like(observed, dtype=int)
    for population in np.unique(block_populations):
        use = block_populations == population
        population_called = called[:, use].astype(float)
        population_dosage = y_filled[:, use]
        population_present = cluster_present[use]
        population_n_called = population_called.sum(axis=1)
        population_with = population_called @ population_present
        population_without = population_n_called[:, None] - population_with
        population_dosage_with = population_dosage @ population_present
        population_dosage_without = (
            population_dosage.sum(axis=1)[:, None] - population_dosage_with
        )
        informative = (
            (population_n_called[:, None] >= min_population_samples)
            & (population_with > 0)
            & (population_without > 0)
        )
        population_effect = np.divide(
            population_dosage_with, population_with,
            out=np.zeros_like(population_with), where=population_with > 0,
        ) - np.divide(
            population_dosage_without, population_without,
            out=np.zeros_like(population_without), where=population_without > 0,
        )
        informative_populations += informative
        consistent_populations += (
            informative & (np.sign(population_effect) == np.sign(observed))
        )
    directional_consistency = np.divide(
        consistent_populations, informative_populations,
        out=np.full_like(observed, np.nan), where=informative_populations > 0,
    )
    directional_consistency[~np.isfinite(observed) | (observed == 0)] = np.nan

    variant_indices, cluster_indices = np.where(eligible)
    output = sv.iloc[variant_positions[variant_indices]][METADATA_COLUMNS].reset_index(
        drop=True
    )
    output["haploblock_id"] = block_id
    output["cluster_id"] = cluster_ids[cluster_indices]
    output["n_called"] = n_called[variant_indices]
    output["n_cluster_carriers"] = n_samples_with_cluster[
        variant_indices, cluster_indices
    ].astype(int)
    output["n_samples_with_cluster"] = output["n_cluster_carriers"]
    output["n_sv_carriers_with_cluster"] = n_sv_carriers_with_cluster[
        variant_indices, cluster_indices
    ].astype(int)
    output["n_sv_noncarriers_with_cluster"] = (
        output["n_samples_with_cluster"] - output["n_sv_carriers_with_cluster"]
    )
    output["cluster_haplotype_count"] = cluster_sum[
        variant_indices, cluster_indices
    ].astype(int)
    output["carrier_rate_with_cluster"] = carrier_rate_with[
        variant_indices, cluster_indices
    ]
    output["carrier_rate_without_cluster"] = carrier_rate_without[
        variant_indices, cluster_indices
    ]
    output["carrier_rate_difference"] = rate_difference[
        variant_indices, cluster_indices
    ]
    output["association_direction"] = direction[variant_indices, cluster_indices]
    output["population_adjusted_r"] = observed[variant_indices, cluster_indices]
    output["p_value"] = p_values[variant_indices, cluster_indices]
    output["permutations_used"] = permutations_used[variant_indices]
    output["informative_populations"] = informative_populations[
        variant_indices, cluster_indices
    ]
    output["directional_consistency"] = directional_consistency[
        variant_indices, cluster_indices
    ]
    return output


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
    threads: int = 1,
    initial_permutations: int | None = None,
    initial_p_threshold: float = 0.1,
    min_refinement_abs_r: float | None = None,
) -> pd.DataFrame:
    samples = [
        sample for sample in metadata["sample_id"].astype(str)
        if sample in sv.columns[len(METADATA_COLUMNS):]
    ]
    sample_metadata = metadata.copy()
    sample_metadata["sample_id"] = sample_metadata["sample_id"].astype(str)
    population_by_sample = sample_metadata.set_index("sample_id")["population"].astype(str).to_dict()
    variant_key = ["sv_record_id"]
    dosage_matrix = dosage_matrix_from_genotypes(sv, samples)
    record_positions = {
        record_id: position for position, record_id in enumerate(sv["sv_record_id"])
    }
    sample_positions = {sample: position for position, sample in enumerate(samples)}
    memberships = memberships[memberships["sample_id"].astype(str).isin(samples)].copy()
    memberships["sample_id"] = memberships["sample_id"].astype(str)
    memberships_by_block = {
        block_id: block_memberships
        for block_id, block_memberships in memberships.groupby("haploblock_id", sort=False)
    }
    block_groups = list(sv_blocks.drop_duplicates(
        [*variant_key, "haploblock_id"]
    ).groupby("haploblock_id", sort=False))
    if not block_groups:
        return pd.DataFrame(columns=ASSOCIATION_COLUMNS)

    def test_block(index_and_group: tuple[int, tuple[str, pd.DataFrame]]) -> list[dict]:
        index, (block_id, block_variants) = index_and_group
        return association_rows_for_block(
            block_id, block_variants,
            memberships_by_block.get(block_id, memberships.iloc[0:0]),
            sv, dosage_matrix, record_positions, sample_positions, population_by_sample,
            permutations, seed + index, min_cluster_haplotypes, min_population_samples,
            refinement_permutations, refinement_p_threshold,
            initial_permutations, initial_p_threshold, min_refinement_abs_r,
        )

    with ThreadPoolExecutor(max_workers=min(threads, len(block_groups))) as executor:
        block_tables = list(executor.map(test_block, enumerate(block_groups)))
    block_tables = [table for table in block_tables if not table.empty]

    if not block_tables:
        return pd.DataFrame(columns=ASSOCIATION_COLUMNS)
    result = assign_pair_q_values(pd.concat(block_tables, ignore_index=True))
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
    depleted = significant & best["association_direction"].eq("carrier_depleted")
    portable = (
        enriched
        & (best["informative_populations"] >= 2)
        & (best["directional_consistency"] >= 0.75)
    )
    dependent = enriched & (best["informative_populations"] >= 2) & ~portable
    best["association_pattern"] = "no_detected_cluster_signal"
    best.loc[depleted, "association_pattern"] = "cluster_exclusion_signal"
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
        "--initial-permutations", type=int, default=200,
        help="Fast triage permutations run before the full screening set.",
    )
    parser.add_argument("--initial-p-threshold", type=float, default=0.1)
    parser.add_argument(
        "--refinement-permutations", type=int, default=1000000,
        help="Independent permutations for pairs passing the refinement threshold.",
    )
    parser.add_argument("--refinement-p-threshold", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-cluster-haplotypes", type=int, default=6)
    parser.add_argument("--min-population-samples", type=int, default=4)
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument(
        "--min-abs-r", type=float, default=0.3,
        help="Minimum absolute effect for permutation advancement and final classification.",
    )
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
                args.refinement_permutations, args.refinement_p_threshold, args.threads,
                args.initial_permutations, args.initial_p_threshold, args.min_abs_r,
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
