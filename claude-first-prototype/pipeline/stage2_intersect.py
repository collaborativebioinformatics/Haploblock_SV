"""Optional Stage 2: classify cohort SVs relative to haploblock boundaries.

Stage 1 already links SVs to associated haplotype clusters. This script is
only needed for the separate positional question. ``position_class`` records
exact overlap: zero blocks is outside, one is within, and two or more is
boundary_crossing. ``near_boundary`` separately records proximity to an edge.
The script reads the same VCF and the compact per-chromosome haploblock tables
registered in Stage 1's config.yaml.
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

from sv_contract import METADATA_COLUMNS, normalize_chrom, parse_info, parse_length, simplify_sv_id


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage2_intersect")


def is_near_boundary(boundaries: np.ndarray, start: int, end: int, boundary_bp: int) -> bool:
    for position in (start, end):
        index = np.searchsorted(boundaries, position)
        if index < len(boundaries) and boundaries[index] - position < boundary_bp:
            return True
        if index and position - boundaries[index - 1] < boundary_bp:
            return True
    return False


def classify_sv_positions(sv: pd.DataFrame, hb: pd.DataFrame, boundary_bp: int) -> pd.DataFrame:
    position_class = np.full(len(sv), "outside_block", dtype=object)
    haploblock_ids = np.full(len(sv), "", dtype=object)
    near_boundary = np.full(len(sv), False, dtype=bool)

    for chrom, hb_group in hb.groupby("chrom", sort=False):
        hb_group = hb_group.sort_values("start")
        starts = hb_group["start"].to_numpy()
        ends = hb_group["end"].to_numpy()
        ids = hb_group["haploblock_id"].to_numpy()
        if len(starts) > 1 and np.any(ends[:-1] > starts[1:]):
            raise ValueError(f"Overlapping haploblocks on {chrom} are not supported")
        boundaries = np.sort(np.concatenate([starts, ends]))
        sv_mask = (sv["chrom"] == chrom).to_numpy()
        if not sv_mask.any():
            continue

        sv_starts = sv.loc[sv_mask, "start"].to_numpy()
        sv_ends = sv.loc[sv_mask, "end"].to_numpy()
        chrom_classes = np.full(sv_mask.sum(), "outside_block", dtype=object)
        chrom_ids = np.full(sv_mask.sum(), "", dtype=object)
        chrom_near_boundary = np.full(sv_mask.sum(), False, dtype=bool)

        for row_index, (sv_start, sv_end) in enumerate(zip(sv_starts, sv_ends)):
            first = np.searchsorted(ends, sv_start, side="right")
            last = np.searchsorted(starts, sv_end, side="left")
            if first < last:
                chrom_classes[row_index] = "within_block" if last - first == 1 else "boundary_crossing"
                chrom_ids[row_index] = ",".join(ids[first:last])
                chrom_near_boundary[row_index] = last - first > 1
            chrom_near_boundary[row_index] |= is_near_boundary(boundaries, sv_start, sv_end, boundary_bp)

        position_class[sv_mask] = chrom_classes
        haploblock_ids[sv_mask] = chrom_ids
        near_boundary[sv_mask] = chrom_near_boundary

    result = sv.copy()
    result["position_class"] = position_class
    result["haploblock_id"] = haploblock_ids
    result["near_boundary"] = near_boundary
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
                    "sv_record_id": f"{chrom}_record_{record_counts[chrom] + 1}",
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
    return {chrom: pd.DataFrame(chrom_rows, columns=METADATA_COLUMNS) for chrom, chrom_rows in rows.items()}


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
            "near_boundary_records": int(annotated["near_boundary"].sum()),
            "output": str(output_path.resolve()),
        }
        log.info("%s: %s", chrom, counts)

    (args.out_dir / "boundary_qc.json").write_text(json.dumps(qc, indent=2) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
