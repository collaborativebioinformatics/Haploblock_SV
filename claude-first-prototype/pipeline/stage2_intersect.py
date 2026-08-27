"""Optional Stage 2: classify cohort SVs relative to haploblock boundaries.

Stage 1 already links SVs to associated haplotype clusters. This script is
only needed for the separate positional question: whether SVs lie inside,
near, or across haploblock boundaries. It reads the same VCF and the compact
per-chromosome haploblock tables registered in Stage 1's config.yaml.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from match_svs_to_clusters import normalize_chrom, parse_info, parse_length, simplify_sv_id


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage2_intersect")


def classify_sv_positions(sv: pd.DataFrame, hb: pd.DataFrame, boundary_bp: int) -> pd.DataFrame:
    position_class = np.full(len(sv), "outside_block", dtype=object)
    haploblock_ids = np.full(len(sv), "", dtype=object)

    for chrom, hb_group in hb.groupby("chrom", sort=False):
        hb_group = hb_group.sort_values("start")
        starts = hb_group["start"].to_numpy()
        ends = hb_group["end"].to_numpy()
        ids = hb_group["haploblock_id"].to_numpy()
        sv_mask = (sv["chrom"] == chrom).to_numpy()
        if not sv_mask.any():
            continue

        sv_starts = sv.loc[sv_mask, "start"].to_numpy()
        sv_ends = sv.loc[sv_mask, "end"].to_numpy()
        chrom_classes = np.full(sv_mask.sum(), "outside_block", dtype=object)
        chrom_ids = np.full(sv_mask.sum(), "", dtype=object)

        for row_index, (sv_start, sv_end) in enumerate(zip(sv_starts, sv_ends)):
            first = np.searchsorted(ends, sv_start - boundary_bp, side="left")
            last = np.searchsorted(starts, sv_end + boundary_bp, side="right")
            relevant = []
            safely_inside = None
            for block_index in range(first, last):
                block_start = starts[block_index]
                block_end = ends[block_index]
                overlaps = sv_start < block_end and sv_end > block_start
                if overlaps:
                    relevant.append(ids[block_index])
                    if (
                        block_start <= sv_start
                        and sv_end <= block_end
                        and sv_start - block_start >= boundary_bp
                        and block_end - sv_end >= boundary_bp
                    ):
                        safely_inside = ids[block_index]
                else:
                    gap = block_start - sv_end if block_start >= sv_end else sv_start - block_end
                    if 0 <= gap < boundary_bp:
                        relevant.append(ids[block_index])

            if not relevant:
                continue
            chrom_classes[row_index] = (
                "within_block" if len(relevant) == 1 and safely_inside is not None else "boundary_crossing"
            )
            chrom_ids[row_index] = ",".join(relevant)

        position_class[sv_mask] = chrom_classes
        haploblock_ids[sv_mask] = chrom_ids

    result = sv.copy()
    result["position_class"] = position_class
    result["haploblock_id"] = haploblock_ids
    return result


def load_sv_metadata(vcf_path: Path, chroms: list[str], max_id_length: int) -> dict[str, pd.DataFrame]:
    selected = set(chroms)
    rows = {chrom: [] for chrom in chroms}
    record_counts: Counter[str] = Counter()
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = normalize_chrom(fields[0])
            if chrom not in selected:
                continue
            position = int(fields[1])
            start = position - 1
            info = parse_info(fields[7])
            end_value = info.get("END")
            end = int(end_value) if isinstance(end_value, str) and end_value != "." else position
            sv_type = str(info.get("SVTYPE", "MISSING"))
            source_id = fields[2] if fields[2] != "." else f"{chrom}:{position}:{sv_type}:{record_counts[chrom] + 1}"
            rows[chrom].append(
                {
                    "sv_id": simplify_sv_id(source_id, chrom, start, end, sv_type, max_id_length),
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "sv_type": sv_type,
                    "length": parse_length(info, start, end),
                    "filter": fields[6],
                    "imprecise": "IMPRECISE" in info,
                }
            )
            record_counts[chrom] += 1
    return {chrom: pd.DataFrame(chrom_rows) for chrom, chrom_rows in rows.items()}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("stage1_output/config.yaml"))
    parser.add_argument("--out-dir", type=Path, default=Path("stage2_boundary_output"))
    parser.add_argument("--boundary-distance-bp", type=int, default=None)
    parser.add_argument("--max-sv-id-length", type=int, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    chroms = list(config["paths"]["haploblocks"])
    boundary_bp = (
        args.boundary_distance_bp
        if args.boundary_distance_bp is not None
        else config["thresholds"]["boundary_distance_bp"]
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)

    max_id_length = (
        args.max_sv_id_length
        if args.max_sv_id_length is not None
        else config["settings"]["max_sv_id_length"]
    )
    sv_by_chrom = load_sv_metadata(Path(config["paths"]["vcf"]), chroms, max_id_length)
    qc = {"boundary_distance_bp": boundary_bp, "chromosomes": {}}
    for chrom in chroms:
        haploblocks = pd.read_csv(config["paths"]["haploblocks"][chrom], sep="\t")
        annotated = classify_sv_positions(sv_by_chrom[chrom], haploblocks, boundary_bp)
        output_path = args.out_dir / f"boundary_svs.{chrom}.tsv"
        annotated.to_csv(output_path, sep="\t", index=False)
        counts = annotated["position_class"].value_counts().to_dict()
        qc["chromosomes"][chrom] = {
            "sv_records": len(annotated),
            "position_class_counts": counts,
            "output": str(output_path.resolve()),
        }
        log.info("%s: %s", chrom, counts)

    (args.out_dir / "boundary_qc.json").write_text(json.dumps(qc, indent=2) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
