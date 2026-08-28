"""Stage 4: classify cohort SVs by population allele frequency.

Reads every chromosome-specific ``sv_genotypes`` table and the normalized
``sample_metadata`` table registered in Stage 1's config. Population labels
are independent of haploblock cluster labels so Stage 6 can compare them
without circularity.

Outputs:
  - sv_af_classification.tsv: one row per SV and population
  - sv_classification.tsv: one row per SV
  - sv_classification_haploblocks.tsv: one row per SV with haploblock assignments
  - haploblock_population_specific_summary.tsv: summary of population-specific SVs within haploblocks
  - stage4_summary.tsv: overall summary of the classification results
  - population_specific_summary.tsv: summary of population-specific SVs across all populations

"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import matplotlib.pyplot as plt

from sv_contract import METADATA_COLUMNS


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage4_classify_af")


def resolve_path(path: str, config_dir: Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else config_dir / path


def genotype_counts(genotype: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    """Return alternate and called allele counts for one sample across all SVs."""
    if pd.api.types.is_numeric_dtype(genotype):
        dosage = pd.to_numeric(genotype, errors="coerce")
        has_call = dosage.notna()
        return dosage.where(has_call, 0).astype(int).to_numpy(), (has_call.astype(int) * 2).to_numpy()

    genotype = genotype.astype("string").str.strip()
    missing = genotype.str.contains(".", regex=False, na=True)
    alternate = genotype.str.count(r"(?:^|[|/])[1-9][0-9]*(?=$|[|/])")
    called = genotype.str.count(r"(?:^|[|/])[0-9]+(?=$|[|/])")
    return (
        alternate.mask(missing, 0).fillna(0).astype(int).to_numpy(),
        called.mask(missing, 0).fillna(0).astype(int).to_numpy(),
    )


def classify(
    sv: pd.DataFrame,
    metadata: pd.DataFrame,
    af_threshold: float,
    absent_af_threshold: float,
    min_samples_per_population: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return separate per-population AF and per-SV classification tables."""
    sv = sv.reset_index(drop=True)
    population_by_sample = dict(
        zip(metadata["sample_id"].astype(str), metadata["population"].astype(str))
    )
    samples = [sample for sample in sv.columns[len(METADATA_COLUMNS):] if sample in population_by_sample]
    samples_by_population = {
        population: [sample for sample in samples if population_by_sample[sample] == population]
        for population in sorted({population_by_sample[sample] for sample in samples})
    }
    population_af = {}
    population_called_samples = {}
    population_called_alleles = {}
    for population, population_samples in samples_by_population.items():
        alt_sum = np.zeros(len(sv), dtype=int)
        allele_sum = np.zeros(len(sv), dtype=int)
        called_samples = np.zeros(len(sv), dtype=int)
        for sample in population_samples:
            alternate, called = genotype_counts(sv[sample])
            alt_sum += alternate
            allele_sum += called
            called_samples += called > 0
        with np.errstate(invalid="ignore", divide="ignore"):
            population_af[population] = np.where(allele_sum > 0, alt_sum / allele_sum, np.nan)
        population_called_samples[population] = called_samples
        population_called_alleles[population] = allele_sum

    classes = []
    specific_populations = []
    other_reasons = []
    populations = list(samples_by_population)
    for row_index in range(len(sv)):
        populations_with_data = [
            population
            for population in populations
            if population_called_samples[population][row_index] >= min_samples_per_population
        ]
        populations_above_threshold = [
            population
            for population in populations_with_data
            if population_af[population][row_index] >= af_threshold
        ]
        other_populations = [
            population
            for population in populations_with_data
            if population not in populations_above_threshold
        ]

        if len(populations_with_data) < 2:
            sv_class, specific_population, reason = "other", "", "insufficient_population_data"
        elif len(populations_above_threshold) >= 2:
            sv_class, specific_population, reason = "common", "", ""
        elif len(populations_above_threshold) == 1 and all(
            population_af[population][row_index] < absent_af_threshold
            for population in other_populations
        ):
            sv_class = "population_specific"
            specific_population = populations_above_threshold[0]
            reason = ""
        elif len(populations_above_threshold) == 1:
            sv_class, specific_population, reason = (
                "other", "", "one_population_high_plus_intermediate_elsewhere"
            )
        else:
            sv_class, specific_population, reason = "other", "", "absent_or_rare"

        classes.append(sv_class)
        specific_populations.append(specific_population)
        other_reasons.append(reason)

    classification_rows = []
    population_rows = []
    for row_index, variant in sv.iterrows():
        classification_rows.append(
            {
                **{column: variant[column] for column in METADATA_COLUMNS},
                "sv_class": classes[row_index],
                "specific_to_population": specific_populations[row_index],
                "other_reason": other_reasons[row_index],
            }
        )
        for population, population_samples in samples_by_population.items():
            af = population_af[population][row_index]
            population_rows.append(
                {
                    **{column: variant[column] for column in METADATA_COLUMNS},
                    "population": population,
                    "n_samples": len(population_samples),
                    "n_called": int(population_called_samples[population][row_index]),
                    "called_alleles": int(population_called_alleles[population][row_index]),
                    "pop_has_data": bool(
                        population_called_samples[population][row_index]
                        >= min_samples_per_population
                    ),
                    "af": float(af) if np.isfinite(af) else np.nan,
                }
            )
    return pd.DataFrame(population_rows), pd.DataFrame(classification_rows)

## New Function to load haploblock assignments
def load_haploblock_assignments(
    config: dict,
    config_dir: Path,
) -> pd.DataFrame:
    """Load one row per SV-haploblock overlap from Stage 1 summaries."""
    tables = []

    for chrom, path in config["paths"]["sv_block_summary"].items():
        block_summary = pd.read_csv(
            resolve_path(path, config_dir),
            sep="\t",
        )

        required_columns = {
            "sv_record_id",
            "haploblock_id",
            "block_start",
            "block_end",
        }
        missing = required_columns - set(block_summary.columns)
        if missing:
            raise ValueError(
                f"{chrom} block summary is missing columns: {sorted(missing)}"
            )

        tables.append(
            block_summary[
                [
                    "sv_record_id",
                    "haploblock_id",
                    "block_start",
                    "block_end",
                ]
            ].drop_duplicates()
        )

    if not tables:
        return pd.DataFrame(
            columns=[
                "sv_record_id",
                "haploblock_id",
                "block_start",
                "block_end",
            ]
        )

    return pd.concat(tables, ignore_index=True).drop_duplicates()

def write_informative_summary(
    classifications: pd.DataFrame,
    classification_haploblocks: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Write summary tables for population-specific SVs and haploblocks."""

    data = classification_haploblocks.copy()
    data["has_haploblock"] = data["haploblock_id"].notna()

    # One row per SV for overall classification counts.
    sv_level = (
        data.groupby("sv_record_id", as_index=False)
        .agg(
            sv_class=("sv_class", "first"),
            specific_to_population=("specific_to_population", "first"),
            haploblock_count=("haploblock_id", "nunique"),
            has_haploblock=("has_haploblock", "any"),
        )
    )

    overall_summary = pd.DataFrame(
        [
            {
                "total_unique_svs": sv_level["sv_record_id"].nunique(),
                "common_svs": int((sv_level["sv_class"] == "common").sum()),
                "population_specific_svs": int(
                    (sv_level["sv_class"] == "population_specific").sum()
                ),
                "other_svs": int((sv_level["sv_class"] == "other").sum()),
                "svs_overlapping_haploblocks": int(
                    sv_level["has_haploblock"].sum()
                ),
                "svs_outside_haploblocks": int(
                    (~sv_level["has_haploblock"]).sum()
                ),
            }
        ]
    )
    overall_summary.to_csv(
        output_dir / "stage4_summary.tsv",
        sep="\t",
        index=False,
    )

    # One row per population-specific SV and population.
    population_summary = (
        data[data["sv_class"] == "population_specific"]
        .dropna(subset=["haploblock_id"])
        .groupby("specific_to_population", as_index=False)
        .agg(
            population_specific_svs=("sv_record_id", "nunique"),
            haploblocks=("haploblock_id", "nunique"),
        )
        .rename(
            columns={
                "specific_to_population": "population",
            }
        )
        .sort_values("population")
    )
    population_summary.to_csv(
        output_dir / "population_specific_summary.tsv",
        sep="\t",
        index=False,
    )

    # One row per haploblock. Use nunique because one SV can produce
    # multiple rows when it overlaps multiple blocks.
    block_data = data.dropna(subset=["haploblock_id"])

    block_summary = (
        block_data.groupby("haploblock_id", as_index=False)
        .agg(
            total_svs=("sv_record_id", "nunique"),
            common_svs=(
                "sv_record_id",
                lambda values: values[
                    block_data.loc[values.index, "sv_class"] == "common"
                ].nunique(),
            ),
            population_specific_svs=(
                "sv_record_id",
                lambda values: values[
                    block_data.loc[values.index, "sv_class"]
                    == "population_specific"
                ].nunique(),
            ),
            populations_with_specific_svs=(
                "specific_to_population",
                lambda values: ", ".join(sorted(set(values.dropna()))),
            ),
            block_start=("block_start", "first"),
            block_end=("block_end", "first"),
        )
        .sort_values(
            ["population_specific_svs", "total_svs"],
            ascending=[False, False],
        )
    )

    block_summary.to_csv(
        output_dir / "haploblock_population_specific_summary.tsv",
        sep="\t",
        index=False,
    )

    log.info("Stage 4 summary:")
    for column, value in overall_summary.iloc[0].items():
        log.info("  %s: %s", column, value)
    log.info(
        "Wrote %s",
        output_dir / "haploblock_population_specific_summary.tsv",
    )
def write_stage4_plots(
    classification_haploblocks: pd.DataFrame,
    output_dir: Path,
    chrom: str,
    top_haploblocks: int,
) -> None:
    """Write chromosome-specific plots for SV class and population specificity."""

    chrom_data = classification_haploblocks[
        classification_haploblocks["chrom"] == chrom
    ].copy()

    chrom_data = chrom_data.dropna(subset=["haploblock_id"])

    if chrom_data.empty:
        log.warning("No haploblock-overlapping SVs available for %s plots", chrom)
        return

    plot_dir = output_dir / chrom
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Count each SV once within each haploblock.
    unique_sv_block = chrom_data.drop_duplicates(
        ["sv_record_id", "haploblock_id"]
    )

    block_counts = (
        unique_sv_block.groupby(
            ["haploblock_id", "sv_class"],
            as_index=False,
        )
        .size()
        .pivot(
            index="haploblock_id",
            columns="sv_class",
            values="size",
        )
        .fillna(0)
    )

    for column in ["common", "population_specific", "other"]:
        if column not in block_counts:
            block_counts[column] = 0

    block_counts = block_counts[["common", "population_specific", "other"]]
    block_counts["total_svs"] = block_counts.sum(axis=1)
    block_counts = block_counts.sort_index()

    # Plot 1: class composition across all haploblocks.
    ax = block_counts[
        ["common", "population_specific", "other"]
    ].plot(
        kind="bar",
        stacked=True,
        figsize=(18, 7),
        color=["#4C78A8", "#E45756", "#B9B9B9"],
        width=0.9,
    )

    ax.set_title(f"{chrom}: SV classification by haploblock")
    ax.set_xlabel("Haploblock")
    ax.set_ylabel("Unique SV count")

    # Show only about 30 readable haploblock labels while retaining
    # every haploblock bar in the plot.
    label_step = max(1, len(block_counts) // 30)
    label_positions = list(range(0, len(block_counts), label_step))
    label_values = [
        str(block_counts.index[position])
        for position in label_positions
    ]

    ax.set_xticks(label_positions)
    ax.set_xticklabels(
        label_values,
        rotation=45,
        ha="right",
        fontsize=7,
    )

    ax.legend(title="SV class")
    ax.figure.tight_layout()
    ax.figure.savefig(
        plot_dir / "haploblock_sv_classification.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(ax.figure)
    
    # Select haploblocks with the largest population-specific SV counts.
    ranked_blocks = (
        block_counts.sort_values(
            ["population_specific", "total_svs"],
            ascending=[False, False],
        )
        .head(top_haploblocks)
        .index
    )

    specific_data = unique_sv_block[
        (unique_sv_block["sv_class"] == "population_specific")
        & unique_sv_block["haploblock_id"].isin(ranked_blocks)
    ]

    # Plot 2: population-specific SVs by haploblock and subpopulation.
    population_counts = (
        specific_data.dropna(subset=["specific_to_population"])
        .drop_duplicates(["sv_record_id", "haploblock_id"])
        .groupby(
            ["haploblock_id", "specific_to_population"],
            as_index=False,
        )
        .size()
        .pivot(
            index="haploblock_id",
            columns="specific_to_population",
            values="size",
        )
        .fillna(0)
    )

    if not population_counts.empty:
        population_counts = population_counts.loc[
            [block for block in ranked_blocks if block in population_counts.index]
        ]

        ax = population_counts.plot(
            kind="bar",
            stacked=True,
            figsize=(18, 8),
            colormap="tab20",
            width=0.9,
        )

        ax.set_title(
            f"{chrom}: population-specific SVs in top haploblocks"
        )
        ax.set_xlabel("Haploblock")
        ax.set_ylabel("Unique population-specific SV count")
        ax.tick_params(axis="x", labelrotation=90)
        ax.legend(
            title="Subpopulation",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
        )
        ax.figure.tight_layout()
        ax.figure.savefig(
            plot_dir / "population_specific_by_haploblock.png",
            dpi=200,
            bbox_inches="tight",
        )
        plt.close(ax.figure)

    # Plot 3: population-specific fraction by genomic position.
    fraction_data = block_counts.copy()
    block_coordinates = (
        unique_sv_block[
            ["haploblock_id", "block_start", "block_end"]
        ]
        .drop_duplicates("haploblock_id")
        .set_index("haploblock_id")
    )

    fraction_data = fraction_data.join(block_coordinates)
    fraction_data["population_specific_fraction"] = (
        fraction_data["population_specific"]
        / fraction_data["total_svs"].replace(0, np.nan)
    )
    fraction_data["midpoint_mb"] = (
        (fraction_data["block_start"] + fraction_data["block_end"]) / 2
    ) / 1_000_000

    ax = fraction_data.plot(
        x="midpoint_mb",
        y="population_specific_fraction",
        kind="scatter",
        figsize=(12, 6),
        s=np.maximum(fraction_data["total_svs"] * 2, 20),
        color="#E45756",
        alpha=0.75,
    )

    ax.set_title(
        f"{chrom}: population-specific SV fraction by haploblock position"
    )
    ax.set_xlabel("Haploblock midpoint (Mb)")
    ax.set_ylabel("Population-specific SVs / total SVs")
    ax.set_ylim(0, 1)
    ax.figure.tight_layout()
    ax.figure.savefig(
        plot_dir / "population_specific_fraction.png",
        dpi=200,
    )
    plt.close(ax.figure)

    block_counts.reset_index().to_csv(
        plot_dir / "haploblock_plot_data.tsv",
        sep="\t",
        index=False,
    )

    if not population_counts.empty:
        population_counts.reset_index().to_csv(
            plot_dir / "population_haploblock_plot_data.tsv",
            sep="\t",
            index=False,
        )

    log.info("Wrote Stage 4 plots for %s to %s", chrom, plot_dir)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("stage1_output/config.yaml"))
    parser.add_argument("--out-dir", type=Path, default=Path("stage4_output"))
    parser.add_argument("--af-threshold", type=float, default=0.05)
    parser.add_argument("--absent-af-threshold", type=float, default=0.01)
    parser.add_argument("--min-samples-per-population", type=int, default=2)
    parser.add_argument( "--plots",action="store_true",help="Generate chromosome-specific Stage 4 plots",)
    parser.add_argument("--plot-dir", type=Path, default=None, help="Directory for plots; defaults to <out-dir>/stage4_plots")
    parser.add_argument("--top-haploblocks", type=int,default=20,help="Number of haploblocks shown in ranked plots")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent
    metadata = pd.read_csv(resolve_path(config["paths"]["sample_metadata"], config_dir), sep="\t")

    population_tables = []
    classification_tables = []
    for chrom, path in config["paths"]["sv_genotypes"].items():
        sv = pd.read_csv(resolve_path(path, config_dir), sep="\t")
        by_population, classifications = classify(
            sv,
            metadata,
            args.af_threshold,
            args.absent_af_threshold,
            args.min_samples_per_population,
        )
        population_tables.append(by_population)
        classification_tables.append(classifications)
        log.info("%s: classified %d SVs", chrom, len(sv))

    by_population = pd.concat(population_tables, ignore_index=True)
    classifications = pd.concat(classification_tables, ignore_index=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    by_population_path = args.out_dir / "sv_af_classification.tsv"
    classification_path = args.out_dir / "sv_classification.tsv"
    by_population.to_csv(by_population_path, sep="\t", index=False)
    classifications.to_csv(classification_path, sep="\t", index=False)

    haploblock_assignments = load_haploblock_assignments(
        config,
        config_dir,
    )
    
    classification_haploblocks = classifications.merge(
        haploblock_assignments,
        on="sv_record_id",
        how="left",
    )
    
    classification_haploblocks_path = (
        args.out_dir / "sv_classification_haploblocks.tsv"
    )
    
    classification_haploblocks.to_csv(
        classification_haploblocks_path,
        sep="\t",
        index=False,
    )

    write_informative_summary(
        classifications,
        classification_haploblocks,
        args.out_dir,
    )

    if args.plots:
        plot_dir = args.plot_dir or args.out_dir / "stage4_plots"

        for chrom in config["paths"]["sv_genotypes"]:
            write_stage4_plots(
                classification_haploblocks,
                plot_dir,
                chrom,
                args.top_haploblocks,
            )

    
    
    class_counts = classifications["sv_class"].value_counts().to_dict()
    log.info("Classified %d SVs: %s", len(classifications), class_counts)
    log.info("Other fraction: %.3f", class_counts.get("other", 0) / len(classifications))

    stage4_config = dict(config)
    stage4_config["thresholds"] = dict(config["thresholds"])
    stage4_config["thresholds"].update(
        {
            "af_common_threshold": args.af_threshold,
            "af_absent_threshold": args.absent_af_threshold,
            "min_samples_per_population": args.min_samples_per_population,
        }
    )
    stage4_config["paths"] = dict(config["paths"])
    stage4_config["paths"].update(
        {
            "sv_af_classification": str(by_population_path.resolve()),
            "sv_classification": str(classification_path.resolve()),
            "sv_classification_haploblocks": str(classification_haploblocks_path.resolve()),
        }
    )
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(stage4_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
