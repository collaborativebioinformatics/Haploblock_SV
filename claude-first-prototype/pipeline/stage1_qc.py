"""Stage 1: QC & normalization for the Haploblock_SV pipeline.

Reads Stage 0's sv_calls.tsv and haploblocks.tsv (paths taken from Stage 0's
config.yaml -- see the contract comment block at the top of stage0_ingest.py),
and:
  1. filters SVs to "high-confidence" calls
  2. drops SVs outside a configurable size range
  3. normalizes breakpoint coordinates
  4. deduplicates exact-duplicate calls
  5. validates the haploblock table is sorted and non-overlapping (raises,
     does not silently continue, on either violation)

Confidence filter, note: dbVar never sets VCF QUAL/FILTER (both are always
"." in nstd152) -- there is no PASS/FAIL field to filter on. The closest
thing nstd152 has to a confidence flag is INFO/IMPRECISE (~28% of real
records), so that is what "high-confidence calls only" means here.

Dedup key, note: naively deduping on (chrom, start, end, sv_type) is wrong
for INS -- dbVar represents every insertion as a ~1bp reference interval
regardless of the inserted sequence's length, so two textually "identical"
INS rows are very often two distinct real insertions that merely share a
locus (verified against real nstd152 data: 2439/50601 rows collide on
(chrom, start, end, sv_type) alone, almost all INS). Including `length` in
the key resolves nearly all of that; the rows that still collide even with
`length` included (612 groups in real nstd152 data, 575 of them one
precise + one imprecise) look like the same event reported twice at
different confidence -- so ties are broken by preferring imprecise=False.

Output contract for pipeline/stage2_intersect.py (not yet implemented):
same TSV format and column set as Stage 0's tables (see stage0_ingest.py's
top-of-file contract), plus:
  - sv_calls.tsv: 0 or more rows removed (confidence/size filters, dedup);
    coordinates satisfy start <= end (swapped if reversed on input, see
    normalize_sv_coordinates()); `imprecise`/`length` columns unchanged.
  - haploblocks.tsv: passed through byte-for-byte unchanged in content
    (Stage 1 only validates it, never edits it) -- copied here so Stage 2
    has one self-contained input directory rather than reaching back into
    Stage 0's output.
  - config.yaml: a copy of Stage 0's config.yaml with `paths` repointed at
    this stage's own sv_calls.tsv/haploblocks.tsv/sample_metadata.tsv, so
    Stage 2 (and beyond) can chain via --config exactly like Stage 1 does.
  - qc_report.json: input/output row counts for every filter step above.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage1_qc")


def load_config(config_path: Path) -> dict:
    with open(config_path) as fh:
        config = yaml.safe_load(fh)
    for key in ("paths", "thresholds", "genome_build"):
        if key not in config:
            raise ValueError(f"{config_path} is missing required top-level key '{key}'")
    return config


# --------------------------------------------------------------------------
# SV filters
# --------------------------------------------------------------------------

def filter_imprecise(sv: pd.DataFrame, drop_imprecise: bool) -> tuple[pd.DataFrame, int]:
    if not drop_imprecise:
        return sv, 0
    if "imprecise" not in sv.columns:
        log.warning("sv_calls has no 'imprecise' column; skipping the confidence filter")
        return sv, 0
    keep = ~sv["imprecise"].astype(bool)
    return sv[keep].copy(), int((~keep).sum())


def filter_size(sv: pd.DataFrame, min_length: int, max_length: int) -> tuple[pd.DataFrame, int, int]:
    if "length" not in sv.columns:
        log.warning("sv_calls has no 'length' column; skipping the size filter")
        return sv, 0, 0
    has_length = sv["length"].notna()
    in_range = (sv["length"] >= min_length) & (sv["length"] <= max_length)
    keep = ~has_length | in_range  # rows with unknown length are exempted, not dropped
    n_out_of_range = int((has_length & ~in_range).sum())
    n_exempted = int((~has_length).sum())
    if n_exempted:
        log.info("%d row(s) have no resolvable length and are exempt from the size filter", n_exempted)
    return sv[keep].copy(), n_out_of_range, n_exempted


def normalize_sv_coordinates(sv: pd.DataFrame, reference_fasta: str | None) -> tuple[pd.DataFrame, int]:
    """Enforce start <= end. True left-alignment (shifting an indel to the
    leftmost coordinate equivalent under the reference sequence, standard
    VCF normalization) needs the reference genome; this prototype does not
    ship or fetch one, so --reference-fasta is accepted but left as an
    explicit no-op rather than shipping an unvalidated implementation.
    """
    if reference_fasta:
        log.warning(
            "--reference-fasta was given but left-alignment is not implemented in this prototype "
            "(needs a validated indel-normalization routine, e.g. bcftools norm -f) -- coordinates "
            "are only sanity-checked (start <= end), not left-aligned"
        )
    sv = sv.copy()
    reversed_mask = sv["start"] > sv["end"]
    n_swapped = int(reversed_mask.sum())
    if n_swapped:
        sv.loc[reversed_mask, ["start", "end"]] = sv.loc[reversed_mask, ["end", "start"]].values
        log.warning("%d row(s) had start > end; swapped", n_swapped)
    return sv, n_swapped


def deduplicate_sv_calls(sv: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Drop exact-duplicate calls, keyed on (chrom, start, end, sv_type, length).

    `length` must be part of the key: dbVar gives every INS the same ~1bp
    reference interval regardless of the inserted sequence's size, so two
    distinct real insertions at the same locus would otherwise look like
    duplicates of each other (see module docstring). Within a duplicate
    group, the precise (imprecise=False) record is kept over an imprecise
    one; ties are broken by sv_id for determinism.
    """
    key_cols = [c for c in ["chrom", "start", "end", "sv_type", "length"] if c in sv.columns]
    sort_cols, ascending = list(key_cols), [True] * len(key_cols)
    if "imprecise" in sv.columns:
        sort_cols.append("imprecise")
        ascending.append(True)  # False (precise) sorts before True (imprecise)
    if "sv_id" in sv.columns:
        sort_cols.append("sv_id")
        ascending.append(True)
    sv_sorted = sv.sort_values(sort_cols, ascending=ascending, na_position="last")
    deduped = sv_sorted.drop_duplicates(subset=key_cols, keep="first")
    n_removed = len(sv) - len(deduped)
    return deduped.sort_index(), n_removed


# --------------------------------------------------------------------------
# Haploblock validation
# --------------------------------------------------------------------------

def validate_haploblocks(hb: pd.DataFrame) -> None:
    """Raise ValueError, listing every offending row, if hb is not sorted by
    start and non-overlapping within each chromosome. Never silently
    re-sorts or drops rows -- an unsorted/overlapping haploblock table is a
    data problem upstream that Stage 2's intersection logic would otherwise
    fail on silently.
    """
    problems = []
    for chrom, group in hb.groupby("chrom", sort=False):
        starts = group["start"].tolist()
        if starts != sorted(starts):
            problems.append(f"{chrom}: rows are not sorted by start as given in the file")

        prev_id, prev_end = None, None
        for _, row in group.sort_values("start").iterrows():
            if prev_end is not None and row["start"] < prev_end:
                problems.append(
                    f"{chrom}: haploblock {prev_id!r} (end={prev_end}) overlaps "
                    f"{row['haploblock_id']!r} (start={row['start']})"
                )
            prev_id, prev_end = row["haploblock_id"], row["end"]

    if problems:
        raise ValueError(
            f"Haploblock BED validation failed ({len(problems)} problem(s)):\n"
            + "\n".join(f"  - {p}" for p in problems)
        )


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", default="example_data/config.yaml", help="Stage 0's config.yaml")
    p.add_argument("--out-dir", default="stage1_output", help="Output directory for cleaned tables, qc_report.json, and this stage's own config.yaml")
    p.add_argument("--min-sv-length", type=int, default=None, help="Overrides config.yaml's thresholds.min_sv_length")
    p.add_argument("--max-sv-length", type=int, default=None, help="Overrides config.yaml's thresholds.max_sv_length")
    p.add_argument("--keep-imprecise", action="store_true", help="Keep IMPRECISE-flagged calls even if config.yaml's thresholds.drop_imprecise is true")
    p.add_argument("--reference-fasta", default=None, help="Reference FASTA for left-alignment (not yet implemented -- see normalize_sv_coordinates())")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    config_path = Path(args.config)
    config = load_config(config_path)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    thresholds = config["thresholds"]
    min_sv_length = args.min_sv_length if args.min_sv_length is not None else thresholds["min_sv_length"]
    max_sv_length = args.max_sv_length if args.max_sv_length is not None else thresholds["max_sv_length"]
    drop_imprecise = bool(thresholds.get("drop_imprecise", True)) and not args.keep_imprecise

    sv = pd.read_csv(config["paths"]["sv_calls"], sep="\t")
    hb = pd.read_csv(config["paths"]["haploblocks"], sep="\t")
    report = {"sv_calls": {"input_rows": len(sv)}, "haploblocks": {"input_rows": len(hb)}}

    sv, n_dropped_imprecise = filter_imprecise(sv, drop_imprecise)
    report["sv_calls"]["dropped_imprecise"] = n_dropped_imprecise
    report["sv_calls"]["rows_after_confidence_filter"] = len(sv)
    log.info("Confidence filter: dropped %d imprecise row(s) (drop_imprecise=%s)", n_dropped_imprecise, drop_imprecise)

    sv, n_out_of_range, n_exempted = filter_size(sv, min_sv_length, max_sv_length)
    report["sv_calls"]["dropped_out_of_size_range"] = n_out_of_range
    report["sv_calls"]["exempted_from_size_filter_unknown_length"] = n_exempted
    report["sv_calls"]["rows_after_size_filter"] = len(sv)
    log.info("Size filter [%d, %d]bp: dropped %d row(s), exempted %d with unknown length", min_sv_length, max_sv_length, n_out_of_range, n_exempted)

    sv, n_swapped = normalize_sv_coordinates(sv, args.reference_fasta)
    report["sv_calls"]["coords_swapped"] = n_swapped

    sv, n_deduped = deduplicate_sv_calls(sv)
    report["sv_calls"]["exact_duplicates_removed"] = n_deduped
    report["sv_calls"]["output_rows"] = len(sv)
    log.info("Dedup: removed %d exact-duplicate row(s)", n_deduped)

    validate_haploblocks(hb)
    report["haploblocks"]["output_rows"] = len(hb)
    report["haploblocks"]["validation"] = "passed (sorted and non-overlapping per chromosome)"
    log.info("Haploblock validation passed: %d block(s), sorted and non-overlapping per chromosome", len(hb))

    sv.to_csv(out_dir / "sv_calls.tsv", sep="\t", index=False)
    hb.to_csv(out_dir / "haploblocks.tsv", sep="\t", index=False)
    shutil.copy(config["paths"]["sample_metadata"], out_dir / "sample_metadata.tsv")

    with open(out_dir / "qc_report.json", "w") as fh:
        json.dump(report, fh, indent=2)

    stage1_config = dict(config)
    stage1_config["paths"] = {
        "sv_calls": str((out_dir / "sv_calls.tsv").resolve()),
        "haploblocks": str((out_dir / "haploblocks.tsv").resolve()),
        "sample_metadata": str((out_dir / "sample_metadata.tsv").resolve()),
    }
    with open(out_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(stage1_config, fh, sort_keys=False)

    log.info(
        "Wrote %d/%d SVs and %d haploblocks to %s",
        len(sv), report["sv_calls"]["input_rows"], len(hb), out_dir.resolve(),
    )
    log.info("QC report: %s", out_dir / "qc_report.json")
    log.info("Config written to %s", (out_dir / "config.yaml").resolve())


if __name__ == "__main__":
    main(sys.argv[1:])
