"""Stage 1: infer chromosome-specific SV-to-haploblock-cluster associations.

This is the primary entry point for the prototype. It runs the cluster-aware
preprocessor on the cohort VCF, writes all products into one Stage 1 output
directory, and records their paths in config.yaml.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import yaml

from match_svs_to_clusters import ALL_CHROMS, request_with_retries, run as run_cluster_aware
from sv_contract import canonical_sample_id, normalize_chrom


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage1_cluster_aware")

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VCF = REPO_ROOT / "input" / "1kgp_ont_cohort.postfilter.full.vcf.gz"
DEFAULT_SAMPLE_METADATA_URL = (
    "https://s3.amazonaws.com/1000g-ont/PROCESSED_DATA/"
    "1kGP_LRSC_500_ONT_Metadata.tsv"
)
METADATA_COLUMN_NAMES = {
    "NHGRI_ID": "sample_id",
    "SubPopulation": "population",
    "SuperPopulation": "superpopulation",
}


def parse_chroms(value: str) -> list[str]:
    if value.strip().lower() == "all":
        return ALL_CHROMS
    return [normalize_chrom(chrom.strip()) for chrom in value.split(",") if chrom.strip()]


def normalize_sample_metadata(source: Path, samples_path: Path, destination: Path) -> None:
    with samples_path.open() as handle:
        sample_rows = list(csv.DictReader(handle, delimiter="\t"))
    original_to_canonical = {
        row["original_sample_id"]: row["sample_id"]
        for row in sample_rows
    }

    with source.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [METADATA_COLUMN_NAMES.get(name, name) for name in reader.fieldnames]
        metadata_by_sample = {}
        for source_row in reader:
            row = {
                METADATA_COLUMN_NAMES.get(name, name): value
                for name, value in source_row.items()
            }
            row["sample_id"] = original_to_canonical.get(
                row["sample_id"], canonical_sample_id(row["sample_id"])
            )
            metadata_by_sample[row["sample_id"]] = row

    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sample in sample_rows:
            if sample["sample_id"] in metadata_by_sample:
                writer.writerow(metadata_by_sample[sample["sample_id"]])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument(
        "--sample-metadata",
        type=Path,
        default=None,
        help="Override the default 1000 Genomes ONT sample metadata table",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("stage1_output"))
    parser.add_argument("--chroms", default="all", help="Comma-separated chromosomes or 'all' for chr1-22,X")
    parser.add_argument("--cluster-base-url", default="https://data.haploblocks.org/haploblock_hashes/1000G")
    parser.add_argument("--cluster-root", type=Path, default=None, help="Use local <root>/<chrom> cluster files")
    parser.add_argument("--chrom-workers", type=int, default=4)
    parser.add_argument("--download-workers", type=int, default=4)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--max-sv-id-length", type=int, default=80)
    parser.add_argument("--association-threshold", type=float, default=0.75)
    parser.add_argument("--posterior-threshold", type=float, default=0.75)
    parser.add_argument("--max-iterations", type=int, default=25)
    parser.add_argument("--tolerance", type=float, default=1e-5)
    parser.add_argument("--boundary-distance-bp", type=int, default=5000)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    chroms = parse_chroms(args.chroms)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    metadata_source = args.sample_metadata
    if metadata_source is None:
        metadata_source = args.out_dir / "1kGP_LRSC_500_ONT_Metadata.tsv"
        log.info("Downloading sample metadata from %s", DEFAULT_SAMPLE_METADATA_URL)
        metadata_source.write_bytes(
            request_with_retries(DEFAULT_SAMPLE_METADATA_URL, args.retries)
        )

    log.info("Running cluster-aware preprocessing for %d chromosome(s)", len(chroms))
    run_cluster_aware(
        argparse.Namespace(
            vcf=args.vcf,
            chroms=",".join(chroms),
            cluster_base_url=args.cluster_base_url,
            cluster_root=args.cluster_root,
            out_dir=args.out_dir,
            chrom_workers=args.chrom_workers,
            download_workers=args.download_workers,
            retries=args.retries,
            max_sv_id_length=args.max_sv_id_length,
            association_threshold=args.association_threshold,
            posterior_threshold=args.posterior_threshold,
            max_iterations=args.max_iterations,
            tolerance=args.tolerance,
        )
    )

    normalized_metadata_path = args.out_dir / "sample_metadata.tsv"
    normalize_sample_metadata(
        metadata_source,
        args.out_dir / "samples.tsv",
        normalized_metadata_path,
    )

    config = {
        "genome_build": "GRCh38",
        "data_sources": {
            "vcf": str(args.vcf.resolve()),
            "haploblock_clusters": (
                str(args.cluster_root.resolve()) if args.cluster_root is not None else args.cluster_base_url
            ),
        },
        "thresholds": {
            "association_probability": args.association_threshold,
            "heterozygote_assignment_probability": args.posterior_threshold,
            "boundary_distance_bp": args.boundary_distance_bp,
        },
        "settings": {
            "max_sv_id_length": args.max_sv_id_length,
            "max_em_iterations": args.max_iterations,
            "em_tolerance": args.tolerance,
        },
        "paths": {
            "vcf": str(args.vcf.resolve()),
            "samples": str((args.out_dir / "samples.tsv").resolve()),
            "sv_genotypes": {
                chrom: str((args.out_dir / f"sv_genotypes.{chrom}.tsv").resolve())
                for chrom in chroms
            },
            "sv_to_clusters": {
                chrom: str((args.out_dir / f"sv_to_clusters.{chrom}.tsv").resolve())
                for chrom in chroms
            },
            "haploblocks": {
                chrom: str((args.out_dir / f"haploblocks.{chrom}.tsv").resolve())
                for chrom in chroms
            },
            "cluster_memberships": {
                chrom: str((args.out_dir / f"cluster_memberships.{chrom}.tsv").resolve())
                for chrom in chroms
            },
            "sv_block_summary": {
                chrom: str((args.out_dir / f"sv_block_summary.{chrom}.tsv").resolve())
                for chrom in chroms
            },
            "debug_and_qc": str((args.out_dir / "debug_and_qc").resolve()),
        },
    }
    config["data_sources"]["sample_metadata"] = (
        str(args.sample_metadata.resolve())
        if args.sample_metadata is not None
        else DEFAULT_SAMPLE_METADATA_URL
    )
    config["paths"]["sample_metadata"] = str(normalized_metadata_path.resolve())
    config_path = args.out_dir / "config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    log.info("Stage 1 complete: %s", config_path.resolve())


if __name__ == "__main__":
    main(sys.argv[1:])
