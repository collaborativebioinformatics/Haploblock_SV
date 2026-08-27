"""Stage 0: merge a cohort of single-sample SV VCFs.

The workflow deliberately separates exact and approximate matching:

1. ``bcftools merge`` unions records from every sample into one multi-sample
   VCF, pasting together per-sample columns.
2. ``truvari collapse`` groups records that are probably the same SV despite
   small representation differences (and, per the note below, does ALL of
   the actual SV matching -- including what would naively look like "exact"
   matches).

No SV matching heuristics are implemented in Python; tune truvari's
comparison thresholds explicitly for each dataset via ``--truvari-arg``.

Why bcftools merge is *not* allowed to match records on its own
-----------------------------------------------------------------
A naive ``bcftools merge`` (the previous version of this script) lets
bcftools decide which records represent "the same" variant using its
default multiallelic-joining behaviour (``-m both``). For SNPs/indels that's
fine, because REF/ALT are fully spelled out. For SVs written with symbolic
ALT alleles (``<DEL>``, ``<INS>``, ``<DUP>``, ...) -- which is how almost
every SV caller writes its output -- bcftools keys its matching off
CHROM/POS/REF/ALT alone and does **not** look at INFO/END or INFO/SVLEN.
Two completely different deletions that merely *start* at the same
position get silently fused into a single multiallelic record, and the
smaller of the two SVLEN/END values quietly wins. This is a documented
bcftools/truvari interaction, not a hypothetical:
https://github.com/ACEnglish/truvari/wiki/collapse#symbolic-variants

Concretely (verified against bcftools 1.19 with two <DEL> records that
share POS=147022730 but have END 147593064 vs. 148013144 -- a 570kb vs. a
990kb deletion):

    bcftools merge (default -m both)   -> ONE row, SV2's true 990kb extent
                                           is silently discarded
    bcftools merge -m id (unique IDs)  -> TWO rows, both extents preserved

The fix used here is the one truvari's own docs recommend: give every
record across the whole cohort a fresh, guaranteed-unique ID (encoding the
sample index, position, and alleles) before merging, then merge with
``-m id``. Because no two input records can ever share an ID, bcftools
merge degrades into a pure "union rows + paste sample columns" operation --
it never joins two records together, so its blind spot around symbolic
SVs simply can't fire. Every scrap of SV matching, including what would
have been an "exact" duplicate, is left to truvari collapse's tunable
comparison engine (refdist/pctseq/pctsize/pctovl), which does look at
END/SVLEN correctly.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import shutil
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage0_vcf_merge")


def read_vcf_manifest(path: Path) -> list[Path]:
    """Read one VCF path per line, ignoring blank lines and ``#`` comments.

    This is intentionally just a text file: to run Stage 0 on a different
    cohort in the future, point --vcf-list at a new manifest. Nothing else
    about the cohort (size, sample names, caller, whether inputs are sorted
    or symbolic-ALT) is assumed here -- that robustness lives in the merge
    step itself, see stage_input_vcf() below.
    """
    if not path.is_file():
        raise FileNotFoundError(f"VCF manifest does not exist: {path}")
    paths: list[Path] = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        value = line.split("#", 1)[0].strip()
        if not value:
            continue
        vcf_path = Path(value).expanduser().resolve()
        if not vcf_path.is_file():
            raise FileNotFoundError(f"{path}:{line_number}: VCF does not exist: {vcf_path}")
        paths.append(vcf_path)
    if not paths:
        raise ValueError(f"VCF manifest contains no input VCFs: {path}")
    if len(paths) < 2:
        raise ValueError("At least two single-sample VCFs are required for cohort merging")
    if len(set(paths)) != len(paths):
        raise ValueError("VCF manifest contains duplicate input paths")
    return paths


def require_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required executable '{name}' was not found on PATH. "
            "Install bcftools and truvari before running Stage 0."
        )


def require_indexed_reference(reference: Path) -> None:
    fai = reference.with_suffix(reference.suffix + ".fai")
    if not reference.is_file():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")
    if not fai.is_file():
        raise FileNotFoundError(
            f"Reference FASTA is not indexed (missing {fai}). Run: samtools faidx {reference}"
        )


def run_command(command: list[str]) -> None:
    log.info("Running: %s", " ".join(command))
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required executable was not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Command failed with exit code {exc.returncode}: {command[0]}") from exc


def capture_command(command: list[str]) -> str:
    log.info("Running: %s", " ".join(command))
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required executable was not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: {command[0]}\n{exc.stderr}"
        ) from exc
    return result.stdout


def get_sample_name(vcf: Path, bcftools: str) -> str:
    names = [n for n in capture_command([bcftools, "query", "-l", str(vcf)]).splitlines() if n]
    if len(names) != 1:
        raise ValueError(f"Expected exactly one sample column in {vcf}, found {len(names)}: {names}")
    return names[0]


def check_sample_names(vcfs: list[Path], bcftools: str) -> list[str]:
    """Fail fast, with a clear message, if two inputs share a sample name.

    bcftools merge will otherwise either error out or silently disambiguate
    (e.g. NAME, NAME.2) depending on --force-samples, which is a confusing
    way to discover that two files in the manifest are the same sample.
    """
    names = [get_sample_name(v, bcftools) for v in vcfs]
    by_name: dict[str, list[Path]] = {}
    for name, path in zip(names, vcfs):
        by_name.setdefault(name, []).append(path)
    duplicates = {name: paths for name, paths in by_name.items() if len(paths) > 1}
    if duplicates:
        detail = "\n".join(
            f"  {name}: {', '.join(str(p) for p in paths)}" for name, paths in duplicates.items()
        )
        raise ValueError(
            "Duplicate sample names across input VCFs:\n"
            f"{detail}\n"
            "Rename the samples first, or pass --allow-duplicate-sample-names to let "
            "bcftools disambiguate them automatically (adds a numeric suffix)."
        )
    return names


def read_contig_lengths(vcf: Path, bcftools: str) -> dict[str, int]:
    """Parse ##contig=<ID=...,length=...> lines out of a VCF header."""
    header = capture_command([bcftools, "view", "-h", str(vcf)])
    lengths: dict[str, int] = {}
    for line in header.splitlines():
        if not line.startswith("##contig="):
            continue
        fields = dict(
            item.split("=", 1) for item in line[len("##contig=<"):-1].split(",") if "=" in item
        )
        if "ID" in fields and "length" in fields:
            try:
                lengths[fields["ID"]] = int(fields["length"])
            except ValueError:
                continue
    return lengths


def require_consistent_contigs(vcfs: list[Path], bcftools: str) -> None:
    """Fail fast if the inputs don't look like they share one reference genome build.

    Cohort merging and truvari collapse both assume every input VCF's POS
    values live on the same coordinate system. Calling the same region
    against different reference builds (e.g. CHM13v2.0 vs. GRCh37 vs.
    GRCh38) -- or even just different chromosome-naming conventions for the
    "same" build (``chr1`` vs. ``1``) -- silently breaks that assumption:
    the same POS means a different physical base pair in each file. This
    tends to surface downstream as a confusing crash (e.g. tabix refusing
    to index truvari's output because chromosome blocks aren't contiguous)
    rather than as an obvious error at the source, so it's checked here
    explicitly, before any staging work is done.
    """
    per_file_contigs = {vcf: read_contig_lengths(vcf, bcftools) for vcf in vcfs}
    shared_names = set.intersection(*(set(c) for c in per_file_contigs.values())) if per_file_contigs else set()

    if not shared_names:
        detail = "\n".join(
            f"  {vcf}: {', '.join(list(contigs)[:5])}{' ...' if len(contigs) > 5 else ''}"
            for vcf, contigs in per_file_contigs.items()
        )
        raise ValueError(
            "Input VCFs share no contig names at all -- they likely use different "
            "chromosome-naming conventions and/or different reference genome builds "
            "(e.g. 'chr1' vs '1', or calls made against CHM13/GRCh37/GRCh38):\n"
            f"{detail}\n"
            "Merging positions from different reference coordinate systems is not "
            "meaningful: the same POS refers to a different physical base pair in "
            "each build. Re-run with VCFs called against a single, shared reference. "
            "If you are certain this is intentional, pass --force-mismatched-references."
        )

    mismatches: dict[str, dict[int, list[Path]]] = {}
    for name in shared_names:
        by_length: dict[int, list[Path]] = {}
        for vcf, contigs in per_file_contigs.items():
            by_length.setdefault(contigs[name], []).append(vcf)
        if len(by_length) > 1:
            mismatches[name] = by_length

    if mismatches:
        detail_lines = []
        for name, by_length in mismatches.items():
            for length, files in by_length.items():
                detail_lines.append(f"  {name} length={length}: {', '.join(str(f) for f in files)}")
        raise ValueError(
            "Input VCFs declare different lengths for the same contig name -- this "
            "almost always means they were called against different reference genome "
            "builds, even though the chromosome names match:\n"
            + "\n".join(detail_lines)
            + "\nRe-run with VCFs called against a single, shared reference. If you are "
            "certain this is intentional, pass --force-mismatched-references."
        )


def stage_input_vcf(vcf: Path, sample_index: int, staging_dir: Path, bcftools: str, threads: int) -> Path:
    """Sort and re-ID one input VCF so it can merge safely with the rest of the cohort.

    Two things happen here, both required before bcftools ever sees more
    than one sample at a time:

    1. Sort. bcftools merge requires every input to be coordinate-sorted;
       an arbitrary future cohort is not guaranteed to already be.
    2. Re-ID. Every record gets a new ID of the form
       ``s<sample_index>_<CHROM>_<POS>_<REF>_<ALT>``, which is unique
       across the whole cohort by construction (the sample-index prefix
       guarantees no two files can collide, even if two samples happen to
       report an identical SV). This is what lets bcftools merge run with
       ``-m id`` -- see the module docstring for why that matters for
       symbolic SV alleles.
    """
    sorted_vcf = staging_dir / f"sample{sample_index:04d}.sorted.vcf.gz"
    staged_vcf = staging_dir / f"sample{sample_index:04d}.staged.vcf.gz"

    run_command([bcftools, "sort", "-Oz", "-o", str(sorted_vcf), str(vcf)])
    run_command(
        [
            bcftools, "annotate",
            "--threads", str(threads),
            "-I", f"s{sample_index}_%CHROM\\_%POS\\_%REF\\_%ALT",
            "-Oz", "-o", str(staged_vcf),
            str(sorted_vcf),
        ]
    )
    run_command([bcftools, "index", "--tbi", "--threads", str(threads), str(staged_vcf)])
    return staged_vcf


def merge_vcfs(
    vcfs: list[Path],
    out_dir: Path,
    truvari_args: list[str] | None = None,
    reference: Path | None = None,
    threads: int = 0,
    allow_duplicate_sample_names: bool = False,
    force_mismatched_references: bool = False,
    bcftools: str = "bcftools",
    truvari: str = "truvari",
) -> dict[str, Path]:
    """Stage every input, bcftools-merge them, then truvari-collapse the result.

    Returns a dict of the artifacts a caller might want: the staged
    per-sample inputs, the exact (post-merge, pre-collapse) VCF, the final
    collapsed VCF, and the removed-output VCF listing what truvari folded
    into each kept record.
    """
    require_executable(bcftools)
    require_executable(truvari)
    if reference is not None:
        require_indexed_reference(reference)

    out_dir.mkdir(parents=True, exist_ok=True)
    staging_dir = out_dir / "staged"
    staging_dir.mkdir(exist_ok=True)

    if not force_mismatched_references:
        require_consistent_contigs(vcfs, bcftools)  # fail fast: same reference build?
    else:
        log.warning(
            "Skipping the contig/reference-build consistency check (--force-mismatched-references). "
            "POS-based matching across mismatched references is not scientifically meaningful."
        )
    if not allow_duplicate_sample_names:
        check_sample_names(vcfs, bcftools)  # fail fast before doing any real work

    staged = [
        stage_input_vcf(vcf, index, staging_dir, bcftools, threads)
        for index, vcf in enumerate(vcfs)
    ]

    merged = out_dir / "merged.exact.vcf.gz"
    merge_command = [bcftools, "merge", "-m", "id", "--threads", str(threads), "-Oz", "-o", str(merged)]
    if allow_duplicate_sample_names:
        merge_command.append("--force-samples")
    merge_command.extend(str(path) for path in staged)
    run_command(merge_command)
    run_command([bcftools, "index", "--tbi", "--threads", str(threads), str(merged)])

    collapsed_raw = out_dir / "merged.collapsed.unsorted.vcf.gz"
    removed_raw = out_dir / "merged.collapsed_removed.unsorted.vcf.gz"
    collapse_command = [
        truvari, "collapse",
        "-i", str(merged),
        "-o", str(collapsed_raw),
        "-c", str(removed_raw),
    ]
    if reference is not None:
        collapse_command.extend(["-f", str(reference)])
    if truvari_args:
        collapse_command.extend(truvari_args)
    run_command(collapse_command)

    # truvari collapse doesn't guarantee its output stays in strict
    # chrom+pos order in every case (e.g. large fractions of the input
    # passed straight through unanalyzed), so re-sort defensively before
    # indexing rather than let tabix fail deep in the pipeline.
    collapsed = out_dir / "merged.collapsed.vcf.gz"
    removed = out_dir / "merged.collapsed_removed.vcf.gz"
    run_command([bcftools, "sort", "-Oz", "-o", str(collapsed), str(collapsed_raw)])
    run_command([bcftools, "sort", "-Oz", "-o", str(removed), str(removed_raw)])
    run_command([bcftools, "index", "--tbi", str(collapsed)])
    run_command([bcftools, "index", "--tbi", str(removed)])

    for label, path in [("merged (pre-collapse)", merged), ("collapsed (kept)", collapsed), ("removed (collapsed-away)", removed)]:
        count = capture_command([bcftools, "index", "-n", str(path)]).strip()
        log.info("%s: %s records -> %s", label, count, path)

    return {
        "staged_inputs": staged,
        "merged": merged,
        "collapsed": collapsed,
        "removed": removed,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--vcf-list",
        type=Path,
        required=True,
        help="Text file containing one local single-sample VCF(.gz) path per line",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("stage0_vcf_merge_output"),
        help="Directory for staged inputs, exact-merge, and collapsed outputs",
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=None,
        help=(
            "Indexed (samtools faidx) reference FASTA used to call the variants. "
            "Passed to truvari collapse as -f/--reference; recommended for datasets "
            "with symbolic SV alleles (<DEL>, <INS>, ...) so truvari can sequence-"
            "resolve them rather than falling back to a coordinate-only estimate."
        ),
    )
    parser.add_argument(
        "--allow-duplicate-sample-names",
        action="store_true",
        help="Let bcftools merge auto-disambiguate duplicate sample names instead of erroring",
    )
    parser.add_argument(
        "--force-mismatched-references",
        action="store_true",
        help=(
            "Skip the preflight check that input VCFs share contig names/lengths "
            "(i.e. were called against the same reference build). Only use this if "
            "you have a specific reason POS-based matching across builds is valid "
            "for your use case -- normally this indicates a real data problem."
        ),
    )
    parser.add_argument("--threads", type=int, default=0, help="Worker threads for bcftools steps")
    parser.add_argument("--bcftools", default="bcftools", help="bcftools executable")
    parser.add_argument("--truvari", default="truvari", help="truvari executable")
    parser.add_argument(
        "--truvari-arg",
        action="append",
        default=[],
        help="Additional argument passed to truvari collapse; repeat as needed",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        vcfs = read_vcf_manifest(args.vcf_list)
        artifacts = merge_vcfs(
            vcfs,
            args.out_dir,
            truvari_args=args.truvari_arg,
            reference=args.reference,
            threads=args.threads,
            allow_duplicate_sample_names=args.allow_duplicate_sample_names,
            force_mismatched_references=args.force_mismatched_references,
            bcftools=args.bcftools,
            truvari=args.truvari,
        )
        run_manifest = {
            "stage": 0,
            "inputs": [str(path) for path in vcfs],
            "reference": str(args.reference) if args.reference else None,
            "force_mismatched_references": args.force_mismatched_references,
            "allow_duplicate_sample_names": args.allow_duplicate_sample_names,
            "staged_inputs": [str(path) for path in artifacts["staged_inputs"]],
            "bcftools_merge_output": str(artifacts["merged"]),
            "truvari_collapse_output": str(artifacts["collapsed"]),
            "truvari_removed_output": str(artifacts["removed"]),
            "truvari_args": args.truvari_arg,
        }
        (args.out_dir / "run.json").write_text(json.dumps(run_manifest, indent=2) + "\n")
        log.info("Wrote exact merge: %s", artifacts["merged"])
        log.info("Wrote collapsed cohort: %s", artifacts["collapsed"])
        log.info("Wrote collapsed-away records: %s", artifacts["removed"])
        return 0
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        log.error("%s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())