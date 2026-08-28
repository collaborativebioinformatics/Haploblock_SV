"""Stage 7: measure whether haploblock hashes resolve SV carriage.

The primary outputs ask whether each resolved SV is concentrated in one or a
few supported clusters and whether each cluster is homogeneous for SV carriage.
Local-diplotype information gain is retained as a secondary description. A
genome-wide SV PCA is written only as population-structure and batch QC.
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from stage4_classify_af import genotype_counts, resolve_path
from sv_contract import METADATA_COLUMNS


COMMON_GT_DOSAGES = {
    "0|0": 0.0, "0/0": 0.0, "0": 0.0,
    "0|1": 1.0, "1|0": 1.0, "0/1": 1.0, "1/0": 1.0, "1": 1.0,
    "1|1": 2.0, "1/1": 2.0,
    ".": np.nan, ".|.": np.nan, "./.": np.nan,
}


def diplotype_table(memberships: pd.DataFrame) -> pd.DataFrame:
    return (
        memberships.sort_values(["haploblock_id", "sample_id", "haplotype"])
        .groupby(["haploblock_id", "sample_id"])["cluster_id"]
        .agg(lambda clusters: "|".join(sorted(map(str, clusters))))
        .rename("diplotype")
        .reset_index()
    )


def dosage_matrix_from_genotypes(sv: pd.DataFrame, samples: list[str]) -> np.ndarray:
    """Return variant-by-sample dosages, using a fast path for ordinary biallelic GTs."""
    matrix = np.empty((len(sv), len(samples)), dtype=float)
    for sample_index, sample in enumerate(samples):
        genotype = sv[sample]
        if pd.api.types.is_numeric_dtype(genotype):
            matrix[:, sample_index] = pd.to_numeric(genotype, errors="coerce")
            continue
        mapped = genotype.map(COMMON_GT_DOSAGES)
        unknown = mapped.isna() & genotype.notna() & ~genotype.isin((".", ".|.", "./."))
        if unknown.any():
            alternate, called = genotype_counts(genotype)
            matrix[:, sample_index] = np.where(called > 0, alternate, np.nan)
        else:
            matrix[:, sample_index] = mapped.to_numpy(dtype=float, na_value=np.nan)
    return matrix


def dosage_table(sv: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """Convert genotype strings once and retain record metadata beside numeric dosages."""
    return pd.concat(
        [
            sv[METADATA_COLUMNS].reset_index(drop=True),
            pd.DataFrame(dosage_matrix_from_genotypes(sv, samples), columns=samples),
        ],
        axis=1,
    )


def binary_entropy_from_counts(carriers: np.ndarray, totals: np.ndarray) -> np.ndarray:
    """Binary entropy for scalar or array carrier/total counts."""
    carriers = np.asarray(carriers, dtype=float)
    totals = np.asarray(totals, dtype=float)
    probability = np.divide(
        carriers, totals, out=np.zeros_like(carriers, dtype=float), where=totals > 0
    )
    complement = 1 - probability
    result = np.zeros_like(probability, dtype=float)
    positive = probability > 0
    result[positive] -= probability[positive] * np.log2(probability[positive])
    positive = complement > 0
    result[positive] -= complement[positive] * np.log2(complement[positive])
    return result


def carrier_cluster_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    """Summarize Stage 1's supported allele-to-cluster assignments."""
    output_columns = [
        *METADATA_COLUMNS, "haploblock_id", "n_supported_carrier_clusters",
        "n_standard_evidence_carrier_clusters", "top_supported_cluster_id",
        "top_cluster_carrier_evidence_share", "effective_carrier_cluster_count",
        "top_standard_evidence_cluster_id", "top_standard_cluster_carrier_evidence_share",
        "effective_standard_carrier_cluster_count",
    ]
    if assignments.empty:
        return pd.DataFrame(columns=output_columns)

    keys = ["sv_record_id", "haploblock_id"]
    unique = assignments.drop_duplicates([*keys, "cluster_id"]).copy()
    unique["_weight"] = pd.to_numeric(
        unique["expected_alt_haplotypes"], errors="coerce"
    ).fillna(0.0)
    unique["_weight_squared"] = np.square(unique["_weight"])

    def summarize(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
        grouped = frame.groupby(keys, sort=False, observed=True)
        totals = grouped.agg(
            cluster_count=("cluster_id", "size"),
            total_weight=("_weight", "sum"),
            squared_weight_sum=("_weight_squared", "sum"),
        ).reset_index()
        top = (
            frame.sort_values("_weight", ascending=False, kind="stable")
            .drop_duplicates(keys)
            [keys + ["cluster_id", "_weight"]]
        )
        result = totals.merge(top, on=keys, how="left")
        positive = result["total_weight"] > 0
        result[f"top_{prefix}_cluster_id"] = result["cluster_id"]
        result[f"top_{prefix}_cluster_carrier_evidence_share"] = np.where(
            positive, result["_weight"] / result["total_weight"], np.nan
        )
        result[f"effective_{prefix}_carrier_cluster_count"] = np.where(
            positive,
            np.square(result["total_weight"]) / result["squared_weight_sum"],
            np.nan,
        )
        return result[
            keys + [
                "cluster_count", f"top_{prefix}_cluster_id",
                f"top_{prefix}_cluster_carrier_evidence_share",
                f"effective_{prefix}_carrier_cluster_count",
            ]
        ]

    supported = summarize(unique, "supported").rename(columns={
        "cluster_count": "n_supported_carrier_clusters",
        "top_supported_cluster_carrier_evidence_share": "top_cluster_carrier_evidence_share",
        "effective_supported_carrier_cluster_count": "effective_carrier_cluster_count",
    })
    standard = unique[unique["evidence_tier"].eq("standard")] if "evidence_tier" in unique else unique
    if standard.empty:
        standard_summary = pd.DataFrame(columns=keys + [
            "n_standard_evidence_carrier_clusters", "top_standard_evidence_cluster_id",
            "top_standard_cluster_carrier_evidence_share",
            "effective_standard_carrier_cluster_count",
        ])
    else:
        standard_summary = summarize(standard, "standard").rename(columns={
            "cluster_count": "n_standard_evidence_carrier_clusters",
            "top_standard_cluster_id": "top_standard_evidence_cluster_id",
        })

    metadata = unique.drop_duplicates(keys)[METADATA_COLUMNS + ["haploblock_id"]]
    result = metadata.merge(supported, on=keys, how="left").merge(
        standard_summary, on=keys, how="left"
    )
    result["n_standard_evidence_carrier_clusters"] = (
        result["n_standard_evidence_carrier_clusters"].fillna(0).astype(int)
    )
    return result[output_columns]


def cluster_population_summary(
    memberships: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Count the populations in which each hash cluster is observed."""
    sample_populations = metadata[["sample_id", "population"]].copy()
    sample_populations["sample_id"] = sample_populations["sample_id"].astype(str)
    represented = memberships.copy()
    represented["sample_id"] = represented["sample_id"].astype(str)
    represented = represented.merge(sample_populations, on="sample_id", how="inner")
    return (
        represented.groupby(["haploblock_id", "cluster_id"])["population"]
        .agg(
            cluster_population_count="nunique",
            cluster_populations=lambda values: ";".join(sorted(set(map(str, values)))),
        )
        .reset_index()
    )


def cluster_purity_table(
    sv: pd.DataFrame,
    sv_blocks: pd.DataFrame,
    memberships: pd.DataFrame,
    samples: list[str],
    min_cluster_samples: int,
    min_carriers: int,
    min_noncarriers: int,
    dosages: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Measure sample-level SV-carriage purity within each haplotype cluster."""
    output_columns = [
        *METADATA_COLUMNS, "haploblock_id", "cluster_id",
        "n_called_cluster_samples", "n_sv_carriers", "n_sv_noncarriers",
        "carrier_rate_in_cluster", "cluster_purity", "mixed_balance",
        "meets_mixed_count_threshold",
    ]
    if dosages is None:
        dosages = dosage_table(sv, samples)
    record_positions = {
        record_id: position
        for position, record_id in enumerate(dosages["sv_record_id"])
    }
    sample_positions = {sample: position for position, sample in enumerate(samples)}
    dosage_values = dosages[samples].to_numpy(dtype=float)
    memberships = memberships[memberships["sample_id"].astype(str).isin(samples)].copy()
    memberships["sample_id"] = memberships["sample_id"].astype(str)
    tables = []

    for block_id, block_variants in sv_blocks.drop_duplicates(
        ["sv_record_id", "haploblock_id"]
    ).groupby("haploblock_id", sort=False):
        block_memberships = memberships[memberships["haploblock_id"] == block_id]
        if block_memberships.empty:
            continue
        cluster_dosage = block_memberships.assign(cluster_haplotypes=1).pivot_table(
            index="sample_id", columns="cluster_id", values="cluster_haplotypes",
            aggfunc="sum", fill_value=0,
        )
        record_ids = [
            record_id for record_id in block_variants["sv_record_id"]
            if record_id in record_positions
        ]
        if not record_ids:
            continue
        variant_positions = np.array(
            [record_positions[record_id] for record_id in record_ids], dtype=int
        )
        block_sample_positions = np.array(
            [sample_positions[sample] for sample in cluster_dosage.index], dtype=int
        )
        block_dosages = dosage_values[
            np.ix_(variant_positions, block_sample_positions)
        ]
        called = np.isfinite(block_dosages)
        carriers = called & (block_dosages > 0)
        cluster_present = cluster_dosage.to_numpy(dtype=bool)
        called_counts = called.astype(np.int32) @ cluster_present.astype(np.int32)
        carrier_counts = carriers.astype(np.int32) @ cluster_present.astype(np.int32)
        variant_indices, cluster_indices = np.where(called_counts >= min_cluster_samples)
        if not len(variant_indices):
            continue

        selected_called = called_counts[variant_indices, cluster_indices]
        selected_carriers = carrier_counts[variant_indices, cluster_indices]
        selected_noncarriers = selected_called - selected_carriers
        carrier_rates = selected_carriers / selected_called
        output = dosages.iloc[variant_positions[variant_indices]][METADATA_COLUMNS].reset_index(
            drop=True
        )
        output["haploblock_id"] = block_id
        output["cluster_id"] = cluster_dosage.columns.to_numpy()[cluster_indices]
        output["n_called_cluster_samples"] = selected_called
        output["n_sv_carriers"] = selected_carriers
        output["n_sv_noncarriers"] = selected_noncarriers
        output["carrier_rate_in_cluster"] = carrier_rates
        output["cluster_purity"] = np.maximum(carrier_rates, 1 - carrier_rates)
        output["mixed_balance"] = 2 * np.minimum(carrier_rates, 1 - carrier_rates)
        output["meets_mixed_count_threshold"] = (
            (selected_carriers >= min_carriers)
            & (selected_noncarriers >= min_noncarriers)
        )
        tables.append(output)
    return (
        pd.concat(tables, ignore_index=True)[output_columns]
        if tables else pd.DataFrame(columns=output_columns)
    )


def representation_summary(
    carrier_clusters: pd.DataFrame,
    purity: pd.DataFrame,
    information: pd.DataFrame,
    purity_threshold: float,
    cluster_populations: pd.DataFrame | None = None,
    classifications: pd.DataFrame | None = None,
    min_cosmopolitan_populations: int = 3,
) -> pd.DataFrame:
    """Combine the SV-centric concentration and cluster-centric purity views."""
    keys = ["sv_record_id", "haploblock_id"]
    output_columns = [
        *METADATA_COLUMNS, "haploblock_id", "n_supported_carrier_clusters",
        "n_standard_evidence_carrier_clusters", "top_supported_cluster_id",
        "top_cluster_carrier_evidence_share", "effective_carrier_cluster_count",
        "top_standard_evidence_cluster_id", "top_standard_cluster_carrier_evidence_share",
        "effective_standard_carrier_cluster_count",
        "n_tested_clusters", "n_mixed_clusters_meeting_count_threshold",
        "n_mixed_diplotypes_meeting_count_threshold",
        "top_supported_cluster_carrier_rate", "top_standard_cluster_carrier_rate",
        "representation_pattern",
        "sv_class", "specific_to_population",
        "top_standard_cluster_population_count", "top_standard_cluster_populations",
        "population_context_pattern",
    ]
    if carrier_clusters.empty and purity.empty and information.empty:
        return pd.DataFrame(columns=output_columns)

    if purity.empty:
        purity_summary = pd.DataFrame(columns=[
            *keys, "n_tested_clusters", "n_mixed_clusters_meeting_count_threshold"
        ])
    else:
        purity_summary = (
            purity.groupby(keys)
            .agg(
                n_tested_clusters=("cluster_id", "size"),
                n_mixed_clusters_meeting_count_threshold=("meets_mixed_count_threshold", "sum"),
            )
            .reset_index()
        )

    result = carrier_clusters.merge(purity_summary, on=keys, how="outer")
    if not information.empty:
        mixed_diplotypes = information[
            keys + ["n_mixed_diplotypes_meeting_count_threshold"]
        ].drop_duplicates(keys)
        result = result.merge(mixed_diplotypes, on=keys, how="outer")
    else:
        result["n_mixed_diplotypes_meeting_count_threshold"] = 0
    if not carrier_clusters.empty and not purity.empty:
        top_rates = (
            carrier_clusters[keys + ["top_supported_cluster_id"]]
            .merge(
                purity[keys + ["cluster_id", "carrier_rate_in_cluster"]],
                left_on=keys + ["top_supported_cluster_id"],
                right_on=keys + ["cluster_id"], how="left",
            )
            [keys + ["carrier_rate_in_cluster"]]
            .rename(columns={"carrier_rate_in_cluster": "top_supported_cluster_carrier_rate"})
        )
        result = result.merge(top_rates, on=keys, how="left")
        standard_top_rates = (
            carrier_clusters[keys + ["top_standard_evidence_cluster_id"]]
            .merge(
                purity[keys + ["cluster_id", "carrier_rate_in_cluster"]],
                left_on=keys + ["top_standard_evidence_cluster_id"],
                right_on=keys + ["cluster_id"], how="left",
            )
            [keys + ["carrier_rate_in_cluster"]]
            .rename(columns={"carrier_rate_in_cluster": "top_standard_cluster_carrier_rate"})
        )
        result = result.merge(standard_top_rates, on=keys, how="left")
    else:
        result["top_supported_cluster_carrier_rate"] = np.nan
        result["top_standard_cluster_carrier_rate"] = np.nan

    if not purity.empty:
        metadata_columns = [column for column in METADATA_COLUMNS if column not in keys]
        purity_metadata = purity.drop_duplicates(keys)[keys + metadata_columns]
        result = result.merge(purity_metadata, on=keys, how="left", suffixes=("", "_purity"))
        for column in metadata_columns:
            purity_column = result.pop(f"{column}_purity") if f"{column}_purity" in result else None
            if column not in result:
                result[column] = purity_column
            elif purity_column is not None:
                result[column] = result[column].fillna(purity_column)
    if not information.empty:
        metadata_columns = [column for column in METADATA_COLUMNS if column not in keys]
        information_metadata = information.drop_duplicates(keys)[keys + metadata_columns]
        result = result.merge(
            information_metadata, on=keys, how="left", suffixes=("", "_information")
        )
        for column in metadata_columns:
            information_column = (
                result.pop(f"{column}_information")
                if f"{column}_information" in result else None
            )
            if column not in result:
                result[column] = information_column
            elif information_column is not None:
                result[column] = result[column].fillna(information_column)
    count_columns = [
        "n_supported_carrier_clusters", "n_standard_evidence_carrier_clusters",
        "n_tested_clusters", "n_mixed_clusters_meeting_count_threshold",
        "n_mixed_diplotypes_meeting_count_threshold",
    ]
    result[count_columns] = result[count_columns].fillna(0).astype(int)
    multiple = result["n_standard_evidence_carrier_clusters"] > 1
    subdivides = result["n_mixed_diplotypes_meeting_count_threshold"] > 0
    single_tag = (
        result["n_standard_evidence_carrier_clusters"].eq(1)
        & result["top_standard_cluster_carrier_rate"].ge(purity_threshold)
        & ~subdivides
    )
    result["representation_pattern"] = "insufficient_or_partial_evidence"
    result.loc[single_tag, "representation_pattern"] = "hash_tag_candidate"
    result.loc[multiple & ~subdivides, "representation_pattern"] = "multi_cluster_sv_candidate"
    result.loc[~multiple & subdivides, "representation_pattern"] = "hash_subdivision_candidate"
    result.loc[multiple & subdivides, "representation_pattern"] = (
        "multi_cluster_and_subdivision_candidate"
    )
    if classifications is not None:
        classification_columns = [
            "sv_record_id", *[
                column for column in ("sv_class", "specific_to_population")
                if column in classifications
            ]
        ]
        result = result.merge(
            classifications[classification_columns].drop_duplicates("sv_record_id"),
            on="sv_record_id", how="left",
        )
    if "sv_class" not in result:
        result["sv_class"] = pd.NA
    if "specific_to_population" not in result:
        result["specific_to_population"] = pd.NA

    if cluster_populations is not None and not cluster_populations.empty:
        population_fields = cluster_populations.rename(columns={
            "cluster_id": "top_standard_evidence_cluster_id",
            "cluster_population_count": "top_standard_cluster_population_count",
            "cluster_populations": "top_standard_cluster_populations",
        })
        result = result.merge(
            population_fields, on=["haploblock_id", "top_standard_evidence_cluster_id"], how="left"
        )
    else:
        result["top_standard_cluster_population_count"] = np.nan
        result["top_standard_cluster_populations"] = pd.NA
    result["population_context_pattern"] = "not_evaluated"
    evaluated = result["sv_class"].notna()
    result.loc[evaluated, "population_context_pattern"] = "no_population_restriction_pattern"
    shared_cluster = result["top_standard_cluster_population_count"].ge(
        min_cosmopolitan_populations
    )
    population_specific = result["sv_class"].eq("population_specific")
    result.loc[
        evaluated & population_specific & shared_cluster,
        "population_context_pattern",
    ] = "population_enriched_on_shared_cluster_candidate"
    return result[output_columns]


def information_table(
    sv: pd.DataFrame,
    sv_blocks: pd.DataFrame,
    memberships: pd.DataFrame,
    samples: list[str],
    min_diplotype_samples: int,
    min_carriers: int = 1,
    min_noncarriers: int = 1,
    dosages: pd.DataFrame | None = None,
    include_information_gain: bool = True,
) -> pd.DataFrame:
    keys = ["sv_record_id"]
    if dosages is None:
        dosages = dosage_table(sv, samples)
    record_positions = {
        record_id: position
        for position, record_id in enumerate(dosages["sv_record_id"])
    }
    sample_positions = {sample: position for position, sample in enumerate(samples)}
    dosage_values = dosages[samples].to_numpy(dtype=float)
    diplotypes = diplotype_table(memberships)
    tables = []

    for block_id, block_variants in sv_blocks.drop_duplicates(
        [*keys, "haploblock_id"]
    ).groupby("haploblock_id", sort=False):
        block_diplotypes = diplotypes[diplotypes["haploblock_id"] == block_id]
        diplotype_by_sample = block_diplotypes.set_index("sample_id")["diplotype"]
        block_samples = [sample for sample in samples if sample in diplotype_by_sample.index]
        record_ids = [
            record_id for record_id in block_variants["sv_record_id"]
            if record_id in record_positions
        ]
        if not block_samples or not record_ids:
            continue

        diplotype_codes, diplotype_labels = pd.factorize(
            diplotype_by_sample.loc[block_samples], sort=True
        )
        indicator = np.eye(len(diplotype_labels), dtype=np.int8)[diplotype_codes]
        variant_positions = np.array(
            [record_positions[record_id] for record_id in record_ids], dtype=int
        )
        block_sample_positions = np.array(
            [sample_positions[sample] for sample in block_samples], dtype=int
        )
        values = dosage_values[np.ix_(variant_positions, block_sample_positions)]
        called = np.isfinite(values)
        carriers = called & (values > 0)
        called_counts = called.astype(np.int32) @ indicator
        carrier_counts = carriers.astype(np.int32) @ indicator
        eligible = called_counts >= min_diplotype_samples
        eligible_counts = np.where(eligible, called_counts, 0)
        eligible_carriers = np.where(eligible, carrier_counts, 0)
        total_samples = eligible_counts.sum(axis=1)
        total_carriers = eligible_carriers.sum(axis=1)
        informative = (total_samples > 0) & (total_carriers > 0) & (
            total_carriers < total_samples
        )
        if not informative.any():
            continue

        noncarrier_counts = called_counts - carrier_counts
        mixed = eligible & (carrier_counts > 0) & (noncarrier_counts > 0)
        mixed_count_threshold = (
            eligible
            & (carrier_counts >= min_carriers)
            & (noncarrier_counts >= min_noncarriers)
        )
        n_eligible = eligible.sum(axis=1)
        output = dosages.iloc[variant_positions][METADATA_COLUMNS].reset_index(drop=True)
        output["haploblock_id"] = block_id
        output["n_samples"] = total_samples
        output["n_diplotypes"] = n_eligible
        output["carrier_rate"] = np.divide(
            total_carriers, total_samples, out=np.zeros_like(total_carriers, dtype=float),
            where=total_samples > 0,
        )
        output["mixed_diplotype_fraction"] = np.divide(
            mixed.sum(axis=1), n_eligible,
            out=np.zeros_like(total_carriers, dtype=float), where=n_eligible > 0,
        )
        output["n_mixed_diplotypes_meeting_count_threshold"] = (
            mixed_count_threshold.sum(axis=1)
        )
        output["samples_in_mixed_diplotypes"] = np.divide(
            np.where(mixed, called_counts, 0).sum(axis=1), total_samples,
            out=np.zeros_like(total_carriers, dtype=float), where=total_samples > 0,
        )

        if include_information_gain:
            baseline_entropy = binary_entropy_from_counts(total_carriers, total_samples)
            group_entropy = binary_entropy_from_counts(carrier_counts, called_counts)
            conditional_entropy = np.divide(
                np.where(eligible, called_counts * group_entropy, 0).sum(axis=1),
                total_samples, out=np.zeros_like(total_carriers, dtype=float),
                where=total_samples > 0,
            )
            information_gain = baseline_entropy - conditional_entropy
            output["carrier_entropy"] = baseline_entropy
            output["conditional_entropy"] = conditional_entropy
            output["information_gain_bits"] = information_gain
            output["normalized_information_gain"] = np.divide(
                information_gain, baseline_entropy,
                out=np.full_like(information_gain, np.nan), where=baseline_entropy > 0,
            )
        else:
            for column in (
                "carrier_entropy", "conditional_entropy", "information_gain_bits",
                "normalized_information_gain",
            ):
                output[column] = np.nan
        tables.append(output.loc[informative])
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def block_summary(information: pd.DataFrame) -> pd.DataFrame:
    if information.empty:
        return pd.DataFrame(columns=[
            "haploblock_id", "n_informative_svs", "mean_normalized_information_gain",
            "mean_mixed_diplotype_fraction", "mean_samples_in_mixed_diplotypes",
        ])
    return (
        information.groupby("haploblock_id")
        .agg(
            n_informative_svs=("sv_id", "size"),
            mean_normalized_information_gain=("normalized_information_gain", "mean"),
            mean_mixed_diplotype_fraction=("mixed_diplotype_fraction", "mean"),
            mean_samples_in_mixed_diplotypes=("samples_in_mixed_diplotypes", "mean"),
        )
        .reset_index()
        .sort_values("mean_samples_in_mixed_diplotypes", ascending=False)
    )


def pca_tables(
    sv: pd.DataFrame,
    metadata: pd.DataFrame,
    min_call_rate: float,
    min_af: float,
    components: int,
    dosage_matrix: np.ndarray | None = None,
    samples: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if dosage_matrix is None:
        samples = [sample for sample in metadata["sample_id"].astype(str) if sample in sv.columns]
        dosage_matrix = dosage_matrix_from_genotypes(sv, samples).T
    elif samples is None:
        raise ValueError("samples are required with a precomputed dosage matrix")

    called = np.isfinite(dosage_matrix)
    called_counts = called.sum(axis=0)
    call_rates = called.mean(axis=0)
    allele_frequencies = np.divide(
        np.nansum(dosage_matrix, axis=0), 2 * called_counts,
        out=np.zeros(dosage_matrix.shape[1], dtype=float), where=called_counts > 0,
    )
    keep = (
        (call_rates >= min_call_rate)
        & (np.minimum(allele_frequencies, 1 - allele_frequencies) >= min_af)
    )
    if not keep.any():
        return pd.DataFrame({"sample_id": samples}), pd.DataFrame(
            columns=["component", "explained_variance_ratio"]
        )

    matrix = dosage_matrix[:, keep].copy()
    frequencies = allele_frequencies[keep]
    matrix[~np.isfinite(matrix)] = np.broadcast_to(
        2 * frequencies, matrix.shape
    )[~np.isfinite(matrix)]
    matrix = (matrix - 2 * frequencies) / np.sqrt(2 * frequencies * (1 - frequencies))
    gram = matrix @ matrix.T
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.clip(eigenvalues[order], 0, None)
    eigenvectors = eigenvectors[:, order]
    n_components = min(components, len(eigenvalues))
    coordinates = eigenvectors[:, :n_components] * np.sqrt(eigenvalues[:n_components])
    coordinate_table = pd.DataFrame(
        coordinates, columns=[f"PC{index + 1}" for index in range(n_components)]
    )
    coordinate_table.insert(0, "sample_id", samples)
    coordinate_table = coordinate_table.merge(metadata, on="sample_id", how="left")
    variance_table = pd.DataFrame({
        "component": [f"PC{index + 1}" for index in range(n_components)],
        "explained_variance_ratio": eigenvalues[:n_components] / eigenvalues.sum(),
        "n_variants": int(keep.sum()),
    })
    return coordinate_table, variance_table


def plot_pca(coordinates: pd.DataFrame, output: Path) -> None:
    if "PC2" not in coordinates:
        return
    import matplotlib.pyplot as plt

    color_fields = [field for field in ("population", "superpopulation") if field in coordinates]
    if not color_fields:
        color_fields = [None]
    figure, axes = plt.subplots(1, len(color_fields), figsize=(6 * len(color_fields), 5), squeeze=False)
    for axis, field in zip(axes[0], color_fields):
        if field is None:
            axis.scatter(coordinates["PC1"], coordinates["PC2"], s=12, alpha=0.7)
            axis.set_title("SV PCA")
        else:
            for label, group in coordinates.groupby(field):
                axis.scatter(group["PC1"], group["PC2"], s=12, alpha=0.7, label=label)
            axis.set_title(f"SV PCA by {field}")
            axis.legend(fontsize=7, frameon=False)
        axis.set_xlabel("PC1")
        axis.set_ylabel("PC2")
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def process_chromosome(
    chrom: str,
    genotype_path: str,
    paths: dict,
    config_dir: Path,
    metadata: pd.DataFrame,
    args: argparse.Namespace,
) -> dict[str, object]:
    """Read and analyze one chromosome using one shared dosage conversion."""
    sv = pd.read_csv(resolve_path(genotype_path, config_dir), sep="\t")
    sv_blocks = pd.read_csv(
        resolve_path(paths["sv_block_summary"][chrom], config_dir), sep="\t"
    )
    memberships = pd.read_csv(
        resolve_path(paths["cluster_memberships"][chrom], config_dir), sep="\t"
    )
    assignments = pd.read_csv(
        resolve_path(paths["sv_to_clusters"][chrom], config_dir), sep="\t"
    )
    samples = [sample for sample in metadata["sample_id"].astype(str) if sample in sv.columns]
    dosages = dosage_table(sv, samples)
    information = information_table(
        sv, sv_blocks, memberships, samples, args.min_diplotype_samples,
        args.min_carriers, args.min_noncarriers, dosages,
        include_information_gain=not args.skip_information_gain,
    )
    purity = cluster_purity_table(
        sv, sv_blocks, memberships, samples, args.min_cluster_samples,
        args.min_carriers, args.min_noncarriers, dosages,
    )
    return {
        "information": information,
        "carrier_clusters": carrier_cluster_summary(assignments),
        "cluster_purity": purity,
        "cluster_populations": cluster_population_summary(memberships, metadata),
        "pca_samples": samples,
        "pca_dosages": None if args.skip_pca else dosages[samples].to_numpy(dtype=float).T,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("stage7_output"))
    parser.add_argument("--min-diplotype-samples", type=int, default=4)
    parser.add_argument("--min-cluster-samples", type=int, default=4)
    parser.add_argument("--min-carriers", type=int, default=3)
    parser.add_argument("--min-noncarriers", type=int, default=3)
    parser.add_argument("--purity-threshold", type=float, default=0.9)
    parser.add_argument("--min-cosmopolitan-populations", type=int, default=3)
    parser.add_argument("--min-call-rate", type=float, default=0.8)
    parser.add_argument("--min-af", type=float, default=0.01)
    parser.add_argument("--components", type=int, default=4)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--skip-information-gain", action="store_true")
    parser.add_argument("--skip-pca", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent
    metadata = pd.read_csv(resolve_path(config["paths"]["sample_metadata"], config_dir), sep="\t")
    classifications = None
    if "sv_classification" in config["paths"]:
        classifications = pd.read_csv(
            resolve_path(config["paths"]["sv_classification"], config_dir), sep="\t"
        )
    chromosome_items = list(config["paths"]["sv_genotypes"].items())
    def analyze(item: tuple[str, str]) -> dict[str, object]:
        return process_chromosome(
            item[0], item[1], config["paths"], config_dir, metadata, args
        )
    if args.threads > 1 and len(chromosome_items) > 1:
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            chromosome_results = list(executor.map(analyze, chromosome_items))
    else:
        chromosome_results = [analyze(item) for item in chromosome_items]

    all_information = [result["information"] for result in chromosome_results]
    all_carrier_clusters = [result["carrier_clusters"] for result in chromosome_results]
    all_cluster_purity = [result["cluster_purity"] for result in chromosome_results]
    all_cluster_populations = [result["cluster_populations"] for result in chromosome_results]
    information = pd.concat(all_information, ignore_index=True)
    summary = block_summary(information)
    carrier_clusters = pd.concat(all_carrier_clusters, ignore_index=True)
    cluster_purity = pd.concat(all_cluster_purity, ignore_index=True)
    cluster_populations = pd.concat(all_cluster_populations, ignore_index=True).drop_duplicates(
        ["haploblock_id", "cluster_id"]
    )
    hash_representation = representation_summary(
        carrier_clusters, cluster_purity, information, args.purity_threshold,
        cluster_populations, classifications, args.min_cosmopolitan_populations,
    )
    if args.skip_pca:
        coordinates = pd.DataFrame()
        variance = pd.DataFrame()
    else:
        pca_samples = list(metadata["sample_id"].astype(str))
        sample_index = {sample: index for index, sample in enumerate(pca_samples)}
        aligned_matrices = []
        for result in chromosome_results:
            matrix = result["pca_dosages"]
            aligned = np.full((len(pca_samples), matrix.shape[1]), np.nan)
            rows = [sample_index[sample] for sample in result["pca_samples"]]
            aligned[rows] = matrix
            aligned_matrices.append(aligned)
        coordinates, variance = pca_tables(
            pd.DataFrame(), metadata, args.min_call_rate, args.min_af, args.components,
            np.concatenate(aligned_matrices, axis=1), pca_samples,
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "sv_haploblock_information": args.out_dir / "sv_haploblock_information.tsv",
        "haploblock_information_summary": args.out_dir / "haploblock_information_summary.tsv",
        "sv_carrier_cluster_summary": args.out_dir / "sv_carrier_cluster_summary.tsv",
        "sv_cluster_purity": args.out_dir / "sv_cluster_purity.tsv",
        "sv_hash_representation": args.out_dir / "sv_hash_representation.tsv",
    }
    if not args.skip_pca:
        paths.update({
            "sv_pca_coordinates": args.out_dir / "sv_pca_coordinates.tsv",
            "sv_pca_variance": args.out_dir / "sv_pca_variance.tsv",
            "sv_pca_plot": args.out_dir / "sv_pca.png",
        })
    information.to_csv(paths["sv_haploblock_information"], sep="\t", index=False)
    summary.to_csv(paths["haploblock_information_summary"], sep="\t", index=False)
    carrier_clusters.to_csv(paths["sv_carrier_cluster_summary"], sep="\t", index=False)
    cluster_purity.to_csv(paths["sv_cluster_purity"], sep="\t", index=False)
    hash_representation.to_csv(paths["sv_hash_representation"], sep="\t", index=False)
    if not args.skip_pca:
        coordinates.to_csv(paths["sv_pca_coordinates"], sep="\t", index=False)
        variance.to_csv(paths["sv_pca_variance"], sep="\t", index=False)
        plot_pca(coordinates, paths["sv_pca_plot"])

    output_config = dict(config)
    output_config["paths"] = dict(config["paths"])
    if args.skip_pca:
        for key in ("sv_pca_coordinates", "sv_pca_variance", "sv_pca_plot"):
            output_config["paths"].pop(key, None)
    output_config["paths"].update({key: str(path.resolve()) for key, path in paths.items()})
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(output_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
