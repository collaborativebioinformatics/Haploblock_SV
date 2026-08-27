"""Stage 2: SV x haploblock intersection for the Haploblock_SV pipeline.

Reads Stage 1's sv_calls.tsv and haploblocks.tsv (paths from --config) and
classifies each SV by how many haploblocks on its chromosome it overlaps
(half-open intervals, so an SV overlaps a block when
`sv.start < block.end and sv.end > block.start`):

  within_block      overlaps exactly one block. (An SV that pokes past the
                     first/last block into the telomeric region still
                     overlaps just that one block -- there is no second
                     haploblock out there to cross into -- so it is
                     within_block, not boundary_crossing.)
  boundary_crossing  overlaps two or more blocks -- i.e. the SV interval
                     physically spans a shared block edge (`sv.start` is
                     below a block's end while `sv.end` is above it / above
                     the next block's start). Blocks are contiguous, so
                     overlapping >1 block is the same thing as straddling a
                     boundary line between them.
  outside_block      overlaps no block on its chromosome.

This is a purely structural test: an SV that sits comfortably inside one
block is `within_block` no matter how close its breakpoints are to that
block's edges. There is no proximity/`boundary_distance_bp` threshold any
more -- being 600 bp from an edge but not crossing it is still
`within_block`. (`--boundary-distance-bp` and `thresholds.boundary_distance_bp`
are kept only so old configs/commands don't error; they are ignored here.)

How the matching is done: for each chromosome, every SV is checked against
*every* haploblock on that chromosome with a plain vectorised interval
comparison (numpy broadcasting) -- no bedtools/pybedtools, and no clever
sorted-search window that could drop a match if the block order were ever
off. Haploblocks are still re-validated as sorted and non-overlapping
(Stage 1 enforces it; Stage 2 re-checks defensively since it is
independently runnable), but correctness of the classification no longer
depends on that ordering.

On `outside_block`: real data.haploblocks.org blocks are contiguous *within
the span they cover* (no inter-block gaps), but that span does not reach
the telomeres/centromere -- e.g. chr21's blocks span ~14.2-46.2 Mb. An SV
before the first block or after the last block on its chromosome is
therefore legitimately `outside_block`; it is not a matching bug. Stage 2
logs a breakdown (before-span / after-span / in an inter-block gap / on a
chromosome with no blocks at all) so a genuine gap problem would be
visible rather than hidden among the expected telomeric calls.

Output contract for pipeline/stage4_classify_af.py (not implemented yet):
same TSV/columns as Stage 1's sv_calls.tsv, plus:
  - position_class (str): one of within_block / boundary_crossing / outside_block
  - haploblock_id   (str): comma-joined ids of every block the SV overlaps
                     (>1 id exactly when position_class is boundary_crossing).
                     In memory this is "" for an outside_block row, but
                     pandas.read_csv reads an empty TSV field back as NaN
                     (float), not "" -- a reader must handle both, e.g.
                     `.fillna("")` before any `.str` accessor call, or
                     filter on position_class instead.
Per-sample genotype columns may be GT strings ("0|1", from Stage 0's --vcf
path) or 0/1/2 dosage ints (dbVar/synthetic paths). haploblocks.tsv and
sample_metadata.tsv are copied through unchanged, and config.yaml is copied
with `paths` repointed at this stage's own output directory, so Stage 4 can
chain via --config the same way Stage 1/2 do.
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
    """Defensive re-check of Stage 1's invariant: per chromosome, blocks are
    sorted by start and non-overlapping. The overlap-count classification
    below does not strictly need the sort, but a violation means the
    haploblock table is malformed -- raise rather than silently mis-classify
    (e.g. counting one physical region twice as two overlapping blocks).
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


def classify_sv_positions(sv: pd.DataFrame, hb: pd.DataFrame) -> pd.DataFrame:
    """Return sv with `position_class` and `haploblock_id` columns added.

    Per chromosome, every SV is compared against every haploblock on that
    chromosome (numpy broadcasting: an [n_sv, n_block] boolean grid). No
    sorted-search shortcut, so the result cannot depend on block ordering.
    Classification is by overlap count only: 0 blocks -> outside_block,
    1 -> within_block, >=2 -> boundary_crossing (the SV physically spans the
    shared edge(s) between the blocks it overlaps).
    """
    position_class = np.full(len(sv), "outside_block", dtype=object)
    haploblock_id = np.full(len(sv), "", dtype=object)
    outside_reason = np.full(len(sv), "no_blocks_on_chrom", dtype=object)

    for chrom, hb_group in hb.groupby("chrom", sort=False):
        hb_group = hb_group.sort_values("start")
        bs = hb_group["start"].to_numpy()
        be = hb_group["end"].to_numpy()
        bid = hb_group["haploblock_id"].to_numpy()
        _validate_sorted_nonoverlapping(chrom, bs, be)

        sv_mask = (sv["chrom"] == chrom).to_numpy()
        if not sv_mask.any():
            continue
        ss = sv.loc[sv_mask, "start"].to_numpy()[:, None]  # [n_sv, 1]
        se = sv.loc[sv_mask, "end"].to_numpy()[:, None]

        # half-open overlap: SV [ss, se) intersects block [bs, be)
        overlaps = (ss < be[None, :]) & (se > bs[None, :])
        n_overlap = overlaps.sum(axis=1)

        chrom_classes = np.where(
            n_overlap == 0,
            "outside_block",
            np.where(n_overlap == 1, "within_block", "boundary_crossing"),
        ).astype(object)
        chrom_block_ids = np.array(
            [",".join(bid[row]) for row in overlaps], dtype=object
        )

        # why is an outside_block SV outside? (blocks exist on this chromosome)
        span_lo, span_hi = bs.min(), be.max()
        se_flat, ss_flat = se[:, 0], ss[:, 0]
        chrom_reason = np.where(
            se_flat <= span_lo, "before_first_block",
            np.where(ss_flat >= span_hi, "after_last_block", "in_inter_block_gap"),
        ).astype(object)

        position_class[sv_mask] = chrom_classes
        haploblock_id[sv_mask] = chrom_block_ids
        outside_reason[sv_mask] = chrom_reason

    out = sv.copy()
    out["position_class"] = position_class
    out["haploblock_id"] = haploblock_id
    # scratch column: logged as a breakdown in main(), then dropped before writing
    out["_outside_reason"] = outside_reason
    return out


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", default="stage1_output/config.yaml", help="Stage 1's config.yaml")
    p.add_argument("--out-dir", default="stage2_output", help="Output directory for the annotated SV table + this stage's own config.yaml")
    p.add_argument("--boundary-distance-bp", type=int, default=None, help="Deprecated / ignored: classification is now structural (SV overlaps >=2 blocks), not proximity-based. Accepted so old commands don't error.")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    config_path = Path(args.config)
    config = load_config(config_path)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.boundary_distance_bp is not None or "boundary_distance_bp" in config.get("thresholds", {}):
        log.info("Note: boundary_distance_bp is no longer used -- an SV is boundary_crossing iff it overlaps >=2 haploblocks")

    sv = pd.read_csv(config["paths"]["sv_calls"], sep="\t")
    hb = pd.read_csv(config["paths"]["haploblocks"], sep="\t")

    annotated = classify_sv_positions(sv, hb)

    counts = annotated["position_class"].value_counts()
    log.info(
        "Position classification: within_block=%d boundary_crossing=%d outside_block=%d",
        int(counts.get("within_block", 0)),
        int(counts.get("boundary_crossing", 0)),
        int(counts.get("outside_block", 0)),
    )
    n_multi_block = (annotated["haploblock_id"].str.count(",") > 0).sum()
    if n_multi_block:
        log.info("%d SV(s) overlap more than one haploblock (span a shared edge)", n_multi_block)

    outside = annotated.loc[annotated["position_class"] == "outside_block", "_outside_reason"]
    if len(outside):
        reasons = outside.value_counts()
        log.info(
            "outside_block breakdown: before_first_block=%d after_last_block=%d "
            "in_inter_block_gap=%d no_blocks_on_chrom=%d",
            int(reasons.get("before_first_block", 0)),
            int(reasons.get("after_last_block", 0)),
            int(reasons.get("in_inter_block_gap", 0)),
            int(reasons.get("no_blocks_on_chrom", 0)),
        )
        n_gap = int(reasons.get("in_inter_block_gap", 0))
        if n_gap:
            log.warning(
                "%d SV(s) fell in an inter-block gap -- data.haploblocks.org blocks are "
                "expected to be contiguous within their span, so a non-zero count here points "
                "at a real gap in the haploblock table, not just telomeric SVs",
                n_gap,
            )
    annotated = annotated.drop(columns="_outside_reason")

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
