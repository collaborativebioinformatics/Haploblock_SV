"""Stage 2: SV x haploblock intersection for the Haploblock_SV pipeline.

Reads Stage 1's sv_calls.tsv and haploblocks.tsv (paths from --config) and
classifies each SV by position relative to the haploblock(s) on its
chromosome:

  within_block      fully contained in exactly one block, and at least
                     `boundary_distance_bp` (N, from config.yaml) away from
                     both of that block's edges.
  boundary_crossing  overlaps more than one block (crosses a shared edge
                     between two blocks), OR is fully contained in one
                     block but within N bp of an edge, OR does not overlap
                     any block but lies within N bp of one (a near-miss,
                     relevant for e.g. a haploblock table with gaps -- real
                     data.haploblocks.org blocks observed so far are
                     contiguous, so this case mostly only fires at the very
                     start/end of a chromosome's covered span).
  outside_block      no block on the SV's chromosome is within N bp at all.

No bedtools/pybedtools: haploblocks are non-overlapping within a
chromosome (Stage 1 enforces this; Stage 2 re-validates it defensively
since it depends on that invariant and must stay independently runnable),
so per chromosome this reduces to a sorted 1D interval search done here
with numpy.searchsorted, not a general-purpose interval-join library.

Output contract for pipeline/stage3_boundary_enrichment.py and
pipeline/stage4_classify_af.py (neither implemented yet): same TSV/columns
as Stage 1's sv_calls.tsv, plus:
  - position_class (str): one of within_block / boundary_crossing / outside_block
  - haploblock_id   (str): comma-joined ids of every block the SV was
                     matched against (see the boundary_crossing definition
                     above for when this is >1 id). In memory this is ""
                     for an outside_block row, but pandas.read_csv reads
                     an empty TSV field back as NaN (float), not "" -- a
                     reader must handle both, e.g. `.fillna("")` before
                     any `.str` accessor call, or filter on position_class
                     instead of on haploblock_id being empty/null.
haploblocks.tsv and sample_metadata.tsv are copied through unchanged (like
Stage 1 does), and config.yaml is copied with `paths` repointed at this
stage's own output directory, so Stage 3/4 can chain via --config the same
way Stage 1/2 do.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage2_intersect")


def load_config(config_path: Path) -> dict:
    with open(config_path) as fh:
        config = yaml.safe_load(fh)
    for key in ("paths", "thresholds"):
        if key not in config:
            raise ValueError(f"{config_path} is missing required top-level key '{key}'")
    return config


def _validate_sorted_nonoverlapping(chrom: str, starts: np.ndarray, ends: np.ndarray) -> None:
    """Defensive re-check of Stage 1's invariant: the searchsorted-based
    matching below is only correct if, per chromosome, blocks are sorted by
    start and non-overlapping. Raise rather than silently mis-classify.
    """
    if len(starts) < 2:
        return
    if not np.all(starts[1:] >= starts[:-1]):
        raise ValueError(f"{chrom}: haploblocks are not sorted by start -- run Stage 1 first")
    if np.any(starts[1:] < ends[:-1]):
        bad = np.nonzero(starts[1:] < ends[:-1])[0]
        raise ValueError(
            f"{chrom}: haploblocks overlap at index/indices {list(bad)} -- run Stage 1 first"
        )


def classify_sv_positions(sv: pd.DataFrame, hb: pd.DataFrame, boundary_bp: int) -> pd.DataFrame:
    """Return sv with `position_class` and `haploblock_id` columns added."""
    position_class = np.full(len(sv), "outside_block", dtype=object)
    haploblock_id = np.full(len(sv), "", dtype=object)

    for chrom, hb_group in hb.groupby("chrom", sort=False):
        hb_group = hb_group.sort_values("start")
        starts = hb_group["start"].to_numpy()
        ends = hb_group["end"].to_numpy()
        ids = hb_group["haploblock_id"].to_numpy()
        _validate_sorted_nonoverlapping(chrom, starts, ends)

        sv_mask = (sv["chrom"] == chrom).to_numpy()
        if not sv_mask.any():
            continue
        sv_starts = sv.loc[sv_mask, "start"].to_numpy()
        sv_ends = sv.loc[sv_mask, "end"].to_numpy()

        # candidate block index range per SV, widened by N so near-misses aren't missed
        idx_lo = np.searchsorted(ends, sv_starts - boundary_bp, side="left")
        idx_hi = np.searchsorted(starts, sv_ends + boundary_bp, side="right") - 1

        chrom_classes = np.full(sv_mask.sum(), "outside_block", dtype=object)
        chrom_block_ids = np.full(sv_mask.sum(), "", dtype=object)
        for row_i, (s, e, lo, hi) in enumerate(zip(sv_starts, sv_ends, idx_lo, idx_hi)):
            if lo > hi:
                continue  # no candidate block within N bp on this chromosome
            relevant = []
            fully_contained_far_from_edges = None
            for j in range(lo, hi + 1):
                b_start, b_end, b_id = starts[j], ends[j], ids[j]
                overlaps = s < b_end and e > b_start
                if overlaps:
                    relevant.append(b_id)
                    contained = b_start <= s and e <= b_end
                    if contained and (s - b_start) >= boundary_bp and (b_end - e) >= boundary_bp:
                        fully_contained_far_from_edges = b_id
                else:
                    gap = (b_start - e) if b_start >= e else (s - b_end)
                    if 0 <= gap < boundary_bp:
                        relevant.append(b_id)
            if not relevant:
                continue
            if len(relevant) == 1 and fully_contained_far_from_edges is not None:
                chrom_classes[row_i] = "within_block"
            else:
                chrom_classes[row_i] = "boundary_crossing"
            chrom_block_ids[row_i] = ",".join(relevant)

        position_class[sv_mask] = chrom_classes
        haploblock_id[sv_mask] = chrom_block_ids

    out = sv.copy()
    out["position_class"] = position_class
    out["haploblock_id"] = haploblock_id
    return out


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", default="stage1_output/config.yaml", help="Stage 1's config.yaml")
    p.add_argument("--out-dir", default="stage2_output", help="Output directory for the annotated SV table + this stage's own config.yaml")
    p.add_argument("--boundary-distance-bp", type=int, default=None, help="Overrides config.yaml's thresholds.boundary_distance_bp")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    config_path = Path(args.config)
    config = load_config(config_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    boundary_bp = args.boundary_distance_bp if args.boundary_distance_bp is not None else config["thresholds"]["boundary_distance_bp"]

    sv = pd.read_csv(config["paths"]["sv_calls"], sep="\t")
    hb = pd.read_csv(config["paths"]["haploblocks"], sep="\t")

    annotated = classify_sv_positions(sv, hb, boundary_bp)

    counts = annotated["position_class"].value_counts()
    log.info(
        "Position classification (N=%dbp): within_block=%d boundary_crossing=%d outside_block=%d",
        boundary_bp,
        int(counts.get("within_block", 0)),
        int(counts.get("boundary_crossing", 0)),
        int(counts.get("outside_block", 0)),
    )
    n_multi_block = (annotated["haploblock_id"].str.count(",") > 0).sum()
    if n_multi_block:
        log.info("%d SV(s) matched more than one haploblock (span a shared edge)", n_multi_block)

    annotated.to_csv(out_dir / "sv_calls.tsv", sep="\t", index=False)
    hb.to_csv(out_dir / "haploblocks.tsv", sep="\t", index=False)
    shutil.copy(config["paths"]["sample_metadata"], out_dir / "sample_metadata.tsv")

    stage2_config = dict(config)
    stage2_config["paths"] = {
        "sv_calls": str((out_dir / "sv_calls.tsv").resolve()),
        "haploblocks": str((out_dir / "haploblocks.tsv").resolve()),
        "sample_metadata": str((out_dir / "sample_metadata.tsv").resolve()),
    }
    with open(out_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(stage2_config, fh, sort_keys=False)

    log.info("Wrote %d annotated SVs and %d haploblocks to %s", len(annotated), len(hb), out_dir.resolve())
    log.info("Config written to %s", (out_dir / "config.yaml").resolve())


if __name__ == "__main__":
    main(sys.argv[1:])
