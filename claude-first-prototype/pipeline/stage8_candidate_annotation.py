"""Stage 8: annotate prioritized SVs with consequence-aware gene evidence.

The same genomic overlap has different interpretations for deletions,
duplications, inversions, and insertions. This stage keeps those consequence
labels explicit and combines them with Stage 6 evidence in a transparent
candidate score intended for triage, not as a claim of causality.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

from stage4_classify_af import resolve_path


def gtf_attributes(text: str) -> dict[str, str]:
    return {
        key: value
        for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text)
    }


def read_gtf(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "exon"}:
                continue
            attributes = gtf_attributes(fields[8])
            rows.append({
                "chrom": fields[0],
                "feature": fields[2],
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "gene_id": attributes.get("gene_id", ""),
                "gene_name": attributes.get("gene_name", attributes.get("gene_id", "")),
                "gene_biotype": attributes.get(
                    "gene_biotype", attributes.get("gene_type", "")
                ),
            })
    return pd.DataFrame(rows)


def overlaps(features: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    return features[(features["start"] < end) & (features["end"] > start)]


def genes_at_breakpoint(features: pd.DataFrame, position: int) -> set[str]:
    genes = features[
        (features["feature"] == "gene")
        & (features["start"] <= position)
        & (features["end"] > position)
    ]
    return set(genes["gene_name"])


def annotate_variant(variant: pd.Series, features: pd.DataFrame) -> tuple[str, str, str]:
    chrom_features = features[features["chrom"] == variant["chrom"]]
    start, end = int(variant["start"]), int(variant["end"])
    sv_type = str(variant["sv_type"]).upper()
    span = overlaps(chrom_features, start, max(end, start + 1))
    genes = set(span[span["feature"] == "gene"]["gene_name"])
    exon_genes = set(span[span["feature"] == "exon"]["gene_name"])

    if sv_type == "DEL":
        if exon_genes:
            return "exon_loss", ";".join(sorted(exon_genes)), "exon"
        if genes:
            return "intragenic_deletion", ";".join(sorted(genes)), "gene"
    elif sv_type == "DUP":
        full_genes = set(
            chrom_features[
                (chrom_features["feature"] == "gene")
                & (chrom_features["start"] >= start)
                & (chrom_features["end"] <= end)
            ]["gene_name"]
        )
        if full_genes:
            return "complete_gene_duplication", ";".join(sorted(full_genes)), "gene"
        if exon_genes or genes:
            affected = exon_genes or genes
            return "partial_gene_duplication", ";".join(sorted(affected)), "exon_or_gene"
    elif sv_type == "INV":
        breakpoint_genes = genes_at_breakpoint(chrom_features, start) | genes_at_breakpoint(
            chrom_features, max(start, end - 1)
        )
        if breakpoint_genes:
            return "inversion_breakpoint_disruption", ";".join(sorted(breakpoint_genes)), "breakpoint"
        if genes:
            return "gene_within_inversion", ";".join(sorted(genes)), "span_only"
    elif sv_type == "INS":
        point = overlaps(chrom_features, start, start + 1)
        point_exons = set(point[point["feature"] == "exon"]["gene_name"])
        point_genes = set(point[point["feature"] == "gene"]["gene_name"])
        if point_exons:
            return "exonic_insertion", ";".join(sorted(point_exons)), "breakpoint"
        if point_genes:
            return "intragenic_insertion", ";".join(sorted(point_genes)), "breakpoint"
    elif exon_genes or genes:
        affected = exon_genes or genes
        return "gene_overlap", ";".join(sorted(affected)), "exon_or_gene"
    return "intergenic_or_unannotated", "", "none"


def consequence_weight(consequence: str) -> int:
    return {
        "exon_loss": 3,
        "exonic_insertion": 3,
        "inversion_breakpoint_disruption": 3,
        "complete_gene_duplication": 2,
        "partial_gene_duplication": 2,
        "intragenic_deletion": 2,
        "intragenic_insertion": 2,
        "gene_overlap": 1,
        "gene_within_inversion": 1,
    }.get(consequence, 0)


def annotate_candidates(
    candidates: pd.DataFrame,
    features: pd.DataFrame,
    classifications: pd.DataFrame | None = None,
) -> pd.DataFrame:
    result = candidates.copy()
    annotations = [annotate_variant(row, features) for _, row in result.iterrows()]
    result[["consequence", "genes", "overlap_basis"]] = pd.DataFrame(
        annotations, index=result.index
    )
    if classifications is not None:
        keys = ["sv_id", "chrom", "start", "end"]
        classification_columns = keys + [
            column for column in ("sv_class", "specific_to_population")
            if column in classifications
        ]
        result = result.merge(
            classifications[classification_columns].drop_duplicates(keys), on=keys, how="left"
        )

    association_weight = result["association_pattern"].map({
        "portable_cluster_tag": 3,
        "cluster_associated": 2,
        "population_dependent_association": 1,
        "no_detected_cluster_signal": 0,
    }).fillna(0)
    statistical_weight = result["q_value"].map(
        lambda value: min(3.0, -math.log10(max(float(value), 1e-300)))
    )
    imprecise = result["imprecise"].astype(str).str.lower().isin({"true", "1"})
    quality_weight = (result["filter"].astype(str) == "PASS").astype(int) + (~imprecise).astype(int)
    result["candidate_score"] = (
        association_weight
        + statistical_weight
        + result["consequence"].map(consequence_weight)
        + quality_weight
    )
    return result.sort_values(
        ["candidate_score", "q_value"], ascending=[False, True]
    ).reset_index(drop=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=Path("stage8_output"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent
    candidates = pd.read_csv(
        resolve_path(config["paths"]["sv_cluster_summary"], config_dir), sep="\t"
    )
    classifications = None
    if "sv_classification" in config["paths"]:
        classifications = pd.read_csv(
            resolve_path(config["paths"]["sv_classification"], config_dir), sep="\t"
        )
    result = annotate_candidates(candidates, read_gtf(args.gtf), classifications)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.out_dir / "annotated_sv_candidates.tsv"
    result.to_csv(output_path, sep="\t", index=False)

    output_config = dict(config)
    output_config["paths"] = dict(config["paths"])
    output_config["paths"]["annotated_sv_candidates"] = str(output_path.resolve())
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(output_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
