"""Stage 7: measure SV information captured or missed by haploblock hashes.

Local diplotypes are used to calculate how much they reduce uncertainty about
SV carriage. A genome-wide SV PCA is also written as descriptive population-
structure and batch QC, not as the primary biological result.
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


def information_table(
    sv: pd.DataFrame,
    sv_blocks: pd.DataFrame,
    memberships: pd.DataFrame,
    samples: list[str],
    min_diplotype_samples: int,
) -> pd.DataFrame:
    keys = ["sv_record_id"]
    dosages = sv[METADATA_COLUMNS].copy()
    for sample in samples:
        alternate, called = genotype_counts(sv[sample])
        dosages[sample] = np.where(called > 0, alternate, np.nan)
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
            for _, group in frame.groupby("diplotype"):
                conditional_entropy += len(group) / len(frame) * entropy(group["carrier"].to_numpy())
                if group["carrier"].nunique() > 1:
                    mixed_groups += 1
                    mixed_samples += len(group)
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
    parser.add_argument("--min-call-rate", type=float, default=0.8)
    parser.add_argument("--min-af", type=float, default=0.01)
    parser.add_argument("--components", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent
    metadata = pd.read_csv(resolve_path(config["paths"]["sample_metadata"], config_dir), sep="\t")
    all_sv = []
    all_information = []
    for chrom, genotype_path in config["paths"]["sv_genotypes"].items():
        sv = pd.read_csv(resolve_path(genotype_path, config_dir), sep="\t")
        sv_blocks = pd.read_csv(
            resolve_path(config["paths"]["sv_block_summary"][chrom], config_dir), sep="\t"
        )
        memberships = pd.read_csv(
            resolve_path(config["paths"]["cluster_memberships"][chrom], config_dir), sep="\t"
        )
        samples = [sample for sample in metadata["sample_id"].astype(str) if sample in sv.columns]
        all_information.append(
            information_table(sv, sv_blocks, memberships, samples, args.min_diplotype_samples)
        )
        all_sv.append(sv)
    information = pd.concat(all_information, ignore_index=True)
    summary = block_summary(information)
    coordinates, variance = pca_tables(
        pd.concat(all_sv, ignore_index=True), metadata,
        args.min_call_rate, args.min_af, args.components,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "sv_haploblock_information": args.out_dir / "sv_haploblock_information.tsv",
        "haploblock_information_summary": args.out_dir / "haploblock_information_summary.tsv",
        "sv_pca_coordinates": args.out_dir / "sv_pca_coordinates.tsv",
        "sv_pca_variance": args.out_dir / "sv_pca_variance.tsv",
        "sv_pca_plot": args.out_dir / "sv_pca.png",
    }
    information.to_csv(paths["sv_haploblock_information"], sep="\t", index=False)
    summary.to_csv(paths["haploblock_information_summary"], sep="\t", index=False)
    coordinates.to_csv(paths["sv_pca_coordinates"], sep="\t", index=False)
    variance.to_csv(paths["sv_pca_variance"], sep="\t", index=False)
    plot_pca(coordinates, paths["sv_pca_plot"])

    output_config = dict(config)
    output_config["paths"] = dict(config["paths"])
    output_config["paths"].update({key: str(path.resolve()) for key, path in paths.items()})
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(output_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
