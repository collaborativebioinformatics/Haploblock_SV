"""Stage 7: measure whether haploblock hashes resolve SV carriage.

The primary outputs ask whether each resolved SV is concentrated in one or a
few supported clusters and whether each cluster is homogeneous for SV carriage.
Local-diplotype information gain is retained as a secondary description. A
genome-wide SV PCA is written only as population-structure and batch QC.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from stage4_classify_af import genotype_counts, resolve_path
from sv_contract import METADATA_COLUMNS


def entropy(values: np.ndarray) -> float:
    if len(values) == 0:
        return np.nan
    frequencies = np.bincount(values.astype(int), minlength=2) / len(values)
    frequencies = frequencies[frequencies > 0]
    return float(-np.sum(frequencies * np.log2(frequencies)))


def diplotype_table(memberships: pd.DataFrame) -> pd.DataFrame:
    return (
        memberships.sort_values(["haploblock_id", "sample_id", "haplotype"])
        .groupby(["haploblock_id", "sample_id"])["cluster_id"]
        .agg(lambda clusters: "|".join(sorted(map(str, clusters))))
        .rename("diplotype")
        .reset_index()
    )


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

    rows = []
    keys = ["sv_record_id", "haploblock_id"]
    for _, group in assignments.drop_duplicates([*keys, "cluster_id"]).groupby(keys, sort=False):
        weights = pd.to_numeric(group["expected_alt_haplotypes"], errors="coerce").fillna(0.0)
        total = float(weights.sum())
        if total > 0:
            shares = weights / total
            top_index = shares.idxmax()
            top_share = float(shares.loc[top_index])
            effective_count = float(1.0 / np.square(shares).sum())
        else:
            top_index = group.index[0]
            top_share = np.nan
            effective_count = np.nan
        standard = group[group["evidence_tier"].eq("standard")] if "evidence_tier" in group else group
        standard_count = len(standard)
        if standard_count:
            standard_weights = pd.to_numeric(
                standard["expected_alt_haplotypes"], errors="coerce"
            ).fillna(0.0)
            standard_total = float(standard_weights.sum())
            if standard_total > 0:
                standard_shares = standard_weights / standard_total
                standard_top_index = standard_shares.idxmax()
                standard_top_cluster = standard.loc[standard_top_index, "cluster_id"]
                standard_top_share = float(standard_shares.loc[standard_top_index])
                standard_effective_count = float(1.0 / np.square(standard_shares).sum())
            else:
                standard_top_cluster = standard.iloc[0]["cluster_id"]
                standard_top_share = np.nan
                standard_effective_count = np.nan
        else:
            standard_top_cluster = pd.NA
            standard_top_share = np.nan
            standard_effective_count = np.nan
        first = group.iloc[0]
        rows.append({
            **{column: first[column] for column in METADATA_COLUMNS},
            "haploblock_id": first["haploblock_id"],
            "n_supported_carrier_clusters": len(group),
            "n_standard_evidence_carrier_clusters": standard_count,
            "top_supported_cluster_id": group.loc[top_index, "cluster_id"],
            "top_cluster_carrier_evidence_share": top_share,
            "effective_carrier_cluster_count": effective_count,
            "top_standard_evidence_cluster_id": standard_top_cluster,
            "top_standard_cluster_carrier_evidence_share": standard_top_share,
            "effective_standard_carrier_cluster_count": standard_effective_count,
        })
    return pd.DataFrame(rows, columns=output_columns)


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
) -> pd.DataFrame:
    """Measure sample-level SV-carriage purity within each haplotype cluster."""
    output_columns = [
        *METADATA_COLUMNS, "haploblock_id", "cluster_id",
        "n_called_cluster_samples", "n_sv_carriers", "n_sv_noncarriers",
        "carrier_rate_in_cluster", "cluster_purity", "mixed_balance",
        "meets_mixed_count_threshold",
    ]
    dosages = sv[METADATA_COLUMNS].copy()
    dosage_columns = {}
    for sample in samples:
        alternate, called = genotype_counts(sv[sample])
        dosage_columns[sample] = np.where(called > 0, alternate, np.nan)
    dosages = pd.concat([dosages, pd.DataFrame(dosage_columns)], axis=1)
    dosages = dosages.set_index("sv_record_id", drop=False)
    memberships = memberships[memberships["sample_id"].astype(str).isin(samples)].copy()
    memberships["sample_id"] = memberships["sample_id"].astype(str)
    rows = []

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
            if record_id in dosages.index
        ]
        if not record_ids:
            continue
        block_dosages = dosages.loc[record_ids, cluster_dosage.index].to_numpy(dtype=float)
        called = np.isfinite(block_dosages)
        carriers = called & (block_dosages > 0)
        cluster_present = cluster_dosage.to_numpy(dtype=bool)
        called_counts = called.astype(int) @ cluster_present.astype(int)
        carrier_counts = carriers.astype(int) @ cluster_present.astype(int)
        for variant_index, record_id in enumerate(record_ids):
            variant = dosages.loc[record_id]
            for cluster_index, cluster_id in enumerate(cluster_dosage.columns):
                n_called = int(called_counts[variant_index, cluster_index])
                if n_called < min_cluster_samples:
                    continue
                n_carriers = int(carrier_counts[variant_index, cluster_index])
                n_noncarriers = n_called - n_carriers
                carrier_rate = n_carriers / n_called
                rows.append({
                    **{column: variant[column] for column in METADATA_COLUMNS},
                    "haploblock_id": block_id,
                    "cluster_id": cluster_id,
                    "n_called_cluster_samples": n_called,
                    "n_sv_carriers": n_carriers,
                    "n_sv_noncarriers": n_noncarriers,
                    "carrier_rate_in_cluster": carrier_rate,
                    "cluster_purity": max(carrier_rate, 1 - carrier_rate),
                    "mixed_balance": 2 * min(carrier_rate, 1 - carrier_rate),
                    "meets_mixed_count_threshold": (
                        n_carriers >= min_carriers and n_noncarriers >= min_noncarriers
                    ),
                })
    return pd.DataFrame(rows, columns=output_columns)


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
) -> pd.DataFrame:
    keys = ["sv_record_id"]
    dosages = sv[METADATA_COLUMNS].copy()
    dosage_columns = {}
    for sample in samples:
        alternate, called = genotype_counts(sv[sample])
        dosage_columns[sample] = np.where(called > 0, alternate, np.nan)
    dosages = pd.concat([dosages, pd.DataFrame(dosage_columns)], axis=1)
    dosages = dosages.set_index(keys, drop=False)
    diplotypes = diplotype_table(memberships)
    rows = []

    for block_id, block_variants in sv_blocks.drop_duplicates(
        [*keys, "haploblock_id"]
    ).groupby("haploblock_id", sort=False):
        block_diplotypes = diplotypes[diplotypes["haploblock_id"] == block_id]
        diplotype_by_sample = block_diplotypes.set_index("sample_id")["diplotype"]
        for _, block_variant in block_variants.iterrows():
            key = block_variant["sv_record_id"]
            if key not in dosages.index:
                continue
            variant = dosages.loc[key]
            called_samples = [
                sample for sample in samples
                if sample in diplotype_by_sample.index and pd.notna(variant[sample])
            ]
            frame = pd.DataFrame({
                "carrier": [int(variant[sample] > 0) for sample in called_samples],
                "diplotype": [diplotype_by_sample[sample] for sample in called_samples],
            })
            counts = frame["diplotype"].value_counts()
            eligible = counts[counts >= min_diplotype_samples].index
            frame = frame[frame["diplotype"].isin(eligible)]
            if frame.empty or frame["carrier"].nunique() < 2:
                continue
            baseline_entropy = entropy(frame["carrier"].to_numpy())
            conditional_entropy = 0.0
            mixed_samples = 0
            mixed_groups = 0
            mixed_groups_meeting_count_threshold = 0
            for _, group in frame.groupby("diplotype"):
                conditional_entropy += len(group) / len(frame) * entropy(group["carrier"].to_numpy())
                if group["carrier"].nunique() > 1:
                    mixed_groups += 1
                    mixed_samples += len(group)
                    n_carriers = int(group["carrier"].sum())
                    n_noncarriers = len(group) - n_carriers
                    if n_carriers >= min_carriers and n_noncarriers >= min_noncarriers:
                        mixed_groups_meeting_count_threshold += 1
            information_gain = baseline_entropy - conditional_entropy
            rows.append({
                **{column: variant[column] for column in METADATA_COLUMNS},
                "haploblock_id": block_id,
                "n_samples": len(frame),
                "n_diplotypes": len(eligible),
                "carrier_rate": float(frame["carrier"].mean()),
                "carrier_entropy": baseline_entropy,
                "conditional_entropy": conditional_entropy,
                "information_gain_bits": information_gain,
                "normalized_information_gain": information_gain / baseline_entropy,
                "mixed_diplotype_fraction": mixed_groups / len(eligible),
                "n_mixed_diplotypes_meeting_count_threshold": (
                    mixed_groups_meeting_count_threshold
                ),
                "samples_in_mixed_diplotypes": mixed_samples / len(frame),
            })
    return pd.DataFrame(rows)


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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = [sample for sample in metadata["sample_id"].astype(str) if sample in sv.columns]
    dosage_matrix = np.empty((len(samples), len(sv)), dtype=float)
    for sample_index, sample in enumerate(samples):
        alternate, called = genotype_counts(sv[sample])
        dosage_matrix[sample_index] = np.where(called > 0, alternate, np.nan)
    columns = []
    for variant_index in range(len(sv)):
        values = dosage_matrix[:, variant_index].copy()
        called = np.isfinite(values)
        if called.mean() < min_call_rate:
            continue
        af = values[called].sum() / (2 * called.sum())
        if min(af, 1 - af) < min_af:
            continue
        values[~called] = values[called].mean()
        scale = np.sqrt(2 * af * (1 - af))
        columns.append((values - 2 * af) / scale)
    if not columns:
        return pd.DataFrame({"sample_id": samples}), pd.DataFrame(
            columns=["component", "explained_variance_ratio"]
        )
    matrix = np.column_stack(columns)
    u, singular_values, _ = np.linalg.svd(matrix, full_matrices=False)
    n_components = min(components, len(singular_values))
    coordinates = u[:, :n_components] * singular_values[:n_components]
    coordinate_table = pd.DataFrame(
        coordinates, columns=[f"PC{index + 1}" for index in range(n_components)]
    )
    coordinate_table.insert(0, "sample_id", samples)
    coordinate_table = coordinate_table.merge(metadata, on="sample_id", how="left")
    variance = singular_values**2
    variance_table = pd.DataFrame({
        "component": [f"PC{index + 1}" for index in range(n_components)],
        "explained_variance_ratio": variance[:n_components] / variance.sum(),
        "n_variants": len(columns),
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
    all_sv = []
    all_information = []
    all_carrier_clusters = []
    all_cluster_purity = []
    all_cluster_populations = []
    for chrom, genotype_path in config["paths"]["sv_genotypes"].items():
        sv = pd.read_csv(resolve_path(genotype_path, config_dir), sep="\t")
        sv_blocks = pd.read_csv(
            resolve_path(config["paths"]["sv_block_summary"][chrom], config_dir), sep="\t"
        )
        memberships = pd.read_csv(
            resolve_path(config["paths"]["cluster_memberships"][chrom], config_dir), sep="\t"
        )
        assignments = pd.read_csv(
            resolve_path(config["paths"]["sv_to_clusters"][chrom], config_dir), sep="\t"
        )
        samples = [sample for sample in metadata["sample_id"].astype(str) if sample in sv.columns]
        all_information.append(
            information_table(
                sv, sv_blocks, memberships, samples, args.min_diplotype_samples,
                args.min_carriers, args.min_noncarriers,
            )
        )
        all_carrier_clusters.append(carrier_cluster_summary(assignments))
        all_cluster_purity.append(cluster_purity_table(
            sv, sv_blocks, memberships, samples, args.min_cluster_samples,
            args.min_carriers, args.min_noncarriers,
        ))
        all_cluster_populations.append(cluster_population_summary(memberships, metadata))
        all_sv.append(sv)
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
    coordinates, variance = pca_tables(
        pd.concat(all_sv, ignore_index=True), metadata,
        args.min_call_rate, args.min_af, args.components,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "sv_haploblock_information": args.out_dir / "sv_haploblock_information.tsv",
        "haploblock_information_summary": args.out_dir / "haploblock_information_summary.tsv",
        "sv_carrier_cluster_summary": args.out_dir / "sv_carrier_cluster_summary.tsv",
        "sv_cluster_purity": args.out_dir / "sv_cluster_purity.tsv",
        "sv_hash_representation": args.out_dir / "sv_hash_representation.tsv",
        "sv_pca_coordinates": args.out_dir / "sv_pca_coordinates.tsv",
        "sv_pca_variance": args.out_dir / "sv_pca_variance.tsv",
        "sv_pca_plot": args.out_dir / "sv_pca.png",
    }
    information.to_csv(paths["sv_haploblock_information"], sep="\t", index=False)
    summary.to_csv(paths["haploblock_information_summary"], sep="\t", index=False)
    carrier_clusters.to_csv(paths["sv_carrier_cluster_summary"], sep="\t", index=False)
    cluster_purity.to_csv(paths["sv_cluster_purity"], sep="\t", index=False)
    hash_representation.to_csv(paths["sv_hash_representation"], sep="\t", index=False)
    coordinates.to_csv(paths["sv_pca_coordinates"], sep="\t", index=False)
    variance.to_csv(paths["sv_pca_variance"], sep="\t", index=False)
    plot_pca(coordinates, paths["sv_pca_plot"])

    output_config = dict(config)
    output_config["paths"] = dict(config["paths"])
    output_config["paths"].update({key: str(path.resolve()) for key, path in paths.items()})
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(output_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
