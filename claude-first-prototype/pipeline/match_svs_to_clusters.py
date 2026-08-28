"""Read an SV VCF and infer haploblock-cluster associations chromosome by chromosome."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import multiprocessing
import re
import subprocess
import threading
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np
import pysam

from sv_contract import (
    METADATA_COLUMNS,
    canonical_sample_id,
    normalize_chrom,
    parse_info,
    parse_length,
    simplify_sv_id,
)


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("match_svs_to_clusters")

_USE_SYSTEM_CURL = False
_DOWNLOAD_CLIENT_LOCK = threading.Lock()
_DOWNLOAD_SESSION = threading.local()

DEFAULT_BASE_URL = "https://data.haploblocks.org/haploblock_hashes/1000G"
ALL_CHROMS = [f"chr{number}" for number in range(1, 23)] + ["chrX"]
CLUSTER_FILENAME_RE = re.compile(r"^(chr[^_]+)_(\d+)-(\d+)_cluster\.tsv$")
HAPLOTYPE_RE = re.compile(
    r"^(?P<sample>.+)_(?P<chrom>chr[^_]+)_region_(?P<start>\d+)-(?P<end>\d+)_hap(?P<haplotype>[01])$"
)


def split_vcf_by_chromosome(
    vcf_path: Path,
    chroms: list[str],
    out_dir: Path,
    max_sv_id_length: int,
    chrom_workers: int = 1,
) -> dict[str, Path]:
    """Write one downstream genotype table per chromosome, using a Tabix index when available."""
    qc_dir = out_dir / "debug_and_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    selected = set(chroms)
    paths = {chrom: out_dir / f"sv_genotypes.{chrom}.tsv" for chrom in chroms}
    handles: dict[str, object] = {}
    writers: dict[str, csv.writer] = {}
    record_counts: Counter[str] = Counter()
    shortened_counts: Counter[str] = Counter()
    filter_counts = {chrom: Counter() for chrom in chroms}
    genotype_counts = {chrom: Counter() for chrom in chroms}
    original_samples: list[str] | None = None
    canonical_samples: list[str] | None = None

    try:
        if not Path(f"{vcf_path}.tbi").exists():
            log.info("Creating a Tabix index for %s", vcf_path)
            pysam.tabix_index(str(vcf_path), preset="vcf")
        with pysam.VariantFile(vcf_path) as source:
            original_samples = list(source.header.samples)
            canonical_samples = [canonical_sample_id(sample) for sample in original_samples]
            if len(canonical_samples) != len(set(canonical_samples)):
                duplicates = sorted(
                    sample for sample, count in Counter(canonical_samples).items() if count > 1
                )
                raise ValueError(f"Sample-ID canonicalization creates duplicate columns: {duplicates}")
            for chrom in chroms:
                handle = paths[chrom].open("w", newline="")
                handles[chrom] = handle
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writers[chrom] = writer
                writer.writerow(METADATA_COLUMNS + canonical_samples)
            def extract_chromosome(requested_chrom: str) -> None:
                with pysam.VariantFile(vcf_path) as indexed_source:
                    for record in indexed_source.fetch(requested_chrom):
                        line = str(record)
                        fields = line.rstrip("\n").split("\t")
                        chrom = normalize_chrom(fields[0])
                        pos = int(fields[1])
                        info = parse_info(fields[7])
                        start = pos - 1
                        end_value = info.get("END")
                        end = int(end_value) if isinstance(end_value, str) and end_value != "." else pos
                        sv_type = str(info.get("SVTYPE", "MISSING"))
                        source_sv_id = fields[2] if fields[2] != "." else f"{chrom}:{pos}:{sv_type}:{record_counts[chrom] + 1}"
                        sv_id = simplify_sv_id(source_sv_id, chrom, start, end, sv_type, max_sv_id_length)
                        shortened_counts[chrom] += sv_id != source_sv_id

                        format_keys = fields[8].split(":")
                        gt_index = format_keys.index("GT")
                        genotypes = []
                        for sample_field in fields[9:]:
                            genotype = sample_field.split(":")[gt_index]
                            genotypes.append(genotype)
                            genotype_counts[chrom][genotype] += 1
                        writers[chrom].writerow(
                            [
                                f"{chrom}_record_{record_counts[chrom] + 1}",
                                sv_id, chrom, start, end, sv_type, parse_length(info, start, end),
                                fields[6], "IMPRECISE" in info,
                            ]
                            + genotypes
                        )
                        record_counts[chrom] += 1
                        filter_counts[chrom][fields[6]] += 1

            with ThreadPoolExecutor(max_workers=min(chrom_workers, len(chroms))) as executor:
                futures = [executor.submit(extract_chromosome, chrom) for chrom in chroms]
                for future in futures:
                    future.result()
    except (OSError, ValueError) as error:
        log.warning("Indexed VCF access unavailable (%s); streaming the VCF instead", error)
        for handle in handles.values():
            handle.close()
        handles.clear()
        writers.clear()
        record_counts.clear()
        shortened_counts.clear()
        filter_counts = {chrom: Counter() for chrom in chroms}
        genotype_counts = {chrom: Counter() for chrom in chroms}
        with gzip.open(vcf_path, "rt") as source:
            for line in source:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    original_samples = line.rstrip("\n").split("\t")[9:]
                    canonical_samples = [canonical_sample_id(sample) for sample in original_samples]
                    if len(canonical_samples) != len(set(canonical_samples)):
                        duplicates = sorted(
                            sample for sample, count in Counter(canonical_samples).items() if count > 1
                        )
                        raise ValueError(f"Sample-ID canonicalization creates duplicate columns: {duplicates}")
                    for chrom in chroms:
                        handle = paths[chrom].open("w", newline="")
                        handles[chrom] = handle
                        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                        writers[chrom] = writer
                        writer.writerow(METADATA_COLUMNS + canonical_samples)
                    continue
                if line.startswith("#"):
                    continue
                if canonical_samples is None:
                    raise ValueError("VCF has no #CHROM header")

                fields = line.rstrip("\n").split("\t")
                chrom = normalize_chrom(fields[0])
                if chrom not in selected:
                    continue
                pos = int(fields[1])
                info = parse_info(fields[7])
                start = pos - 1
                end_value = info.get("END")
                end = int(end_value) if isinstance(end_value, str) and end_value != "." else pos
                sv_type = str(info.get("SVTYPE", "MISSING"))
                source_sv_id = fields[2] if fields[2] != "." else f"{chrom}:{pos}:{sv_type}:{record_counts[chrom] + 1}"
                sv_id = simplify_sv_id(source_sv_id, chrom, start, end, sv_type, max_sv_id_length)
                shortened_counts[chrom] += sv_id != source_sv_id

                format_keys = fields[8].split(":")
                gt_index = format_keys.index("GT")
                genotypes = []
                for sample_field in fields[9:]:
                    genotype = sample_field.split(":")[gt_index]
                    genotypes.append(genotype)
                    genotype_counts[chrom][genotype] += 1
                writers[chrom].writerow(
                    [
                        f"{chrom}_record_{record_counts[chrom] + 1}",
                        sv_id, chrom, start, end, sv_type, parse_length(info, start, end),
                        fields[6], "IMPRECISE" in info,
                    ]
                    + genotypes
                )
                record_counts[chrom] += 1
                filter_counts[chrom][fields[6]] += 1
    finally:
        for handle in handles.values():
            handle.close()

    if original_samples is None or canonical_samples is None:
        raise ValueError("VCF has no #CHROM header")
    with (out_dir / "samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "original_sample_id"])
        writer.writerows(zip(canonical_samples, original_samples))
    for chrom in chroms:
        qc = {
            "vcf": str(vcf_path.resolve()),
            "chrom": chrom,
            "records_written": record_counts[chrom],
            "sv_id_max_length": max_sv_id_length,
            "sv_ids_shortened": shortened_counts[chrom],
            "samples": len(canonical_samples),
            "sample_ids_changed": sum(a != b for a, b in zip(original_samples, canonical_samples)),
            "filter_counts": dict(sorted(filter_counts[chrom].items())),
            "genotype_counts": dict(sorted(genotype_counts[chrom].items())),
        }
        (qc_dir / f"vcf_qc.{chrom}.json").write_text(json.dumps(qc, indent=2) + "\n")
        log.info("Extracted %d records for %s", record_counts[chrom], chrom)
    return paths


def reusable_sv_tables(
    vcf_path: Path,
    chroms: list[str],
    out_dir: Path,
    max_sv_id_length: int,
) -> dict[str, Path] | None:
    """Return a complete interrupted-run VCF split when its QC matches this invocation."""
    qc_dir = out_dir / "debug_and_qc"
    tables = {chrom: out_dir / f"sv_genotypes.{chrom}.tsv" for chrom in chroms}
    if not (out_dir / "samples.tsv").exists():
        return None
    expected_vcf = str(vcf_path.resolve())
    for chrom, table in tables.items():
        qc_path = qc_dir / f"vcf_qc.{chrom}.json"
        if not table.exists() or not qc_path.exists():
            return None
        qc = json.loads(qc_path.read_text())
        if qc.get("vcf") != expected_vcf or qc.get("sv_id_max_length") != max_sv_id_length:
            return None
    log.info("Reusing the completed chromosome-specific VCF split from an interrupted run")
    return tables


@dataclass
class SVRecord:
    metadata: list[str]
    start: int
    end: int
    dosages: np.ndarray


@dataclass
class BlockMembership:
    haploblock_id: str
    chrom: str
    start: int
    end: int
    sample_clusters: dict[str, tuple[str, str]]
    total_cluster_haplotypes: Counter[str]


@dataclass
class PreparedBlock:
    membership: BlockMembership
    sample_indices: np.ndarray
    cluster_ids: list[str]
    cluster0: np.ndarray
    cluster1: np.ndarray
    represented_counts: np.ndarray


def download_with_system_curl(url: str, timeout: int) -> bytes:
    """Download with the OS HTTPS client, which trusts the macOS keychain."""
    result = subprocess.run(
        ["curl", "--fail", "--silent", "--show-error", "--location", "--max-time", str(timeout), url],
        check=True,
        capture_output=True,
    )
    return result.stdout


def request_with_retries(url: str, retries: int, timeout: int = 60) -> bytes:
    import requests

    global _USE_SYSTEM_CURL
    error: Exception | None = None
    for attempt in range(retries):
        try:
            with _DOWNLOAD_CLIENT_LOCK:
                use_system_curl = _USE_SYSTEM_CURL
            if use_system_curl:
                return download_with_system_curl(url, timeout)
            session = getattr(_DOWNLOAD_SESSION, "session", None)
            if session is None:
                session = requests.Session()
                _DOWNLOAD_SESSION.session = session
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            return response.content
        except requests.RequestException as exc:
            with _DOWNLOAD_CLIENT_LOCK:
                if not _USE_SYSTEM_CURL:
                    log.warning(
                        "Python's HTTPS client could not reach the server; using the system HTTPS client instead"
                    )
                    _USE_SYSTEM_CURL = True
            try:
                return download_with_system_curl(url, timeout)
            except subprocess.CalledProcessError as exc:
                error = exc
        except subprocess.CalledProcessError as exc:
            error = exc
        if attempt + 1 < retries:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed to download {url} after {retries} attempts: {error}")


def discover_cluster_filenames(base_url: str, chrom: str, retries: int) -> list[str]:
    directory_url = f"{base_url.rstrip('/')}/{chrom}/clusters/"
    log.info("Discovering cluster files for %s", chrom)
    html = request_with_retries(directory_url, retries).decode()
    pattern = re.compile(rf"{re.escape(chrom)}_\d+-\d+_cluster\.tsv")
    names = sorted(set(pattern.findall(html)), key=parse_cluster_filename)
    if not names:
        raise RuntimeError(f"No cluster files found at {directory_url}")
    return names


def parse_cluster_filename(path_or_name: str | Path) -> tuple[str, int, int]:
    name = Path(path_or_name).name
    match = CLUSTER_FILENAME_RE.fullmatch(name)
    if not match:
        raise ValueError(f"Unexpected cluster filename: {name}")
    return match.group(1), int(match.group(2)), int(match.group(3))


def cluster_cache_directory(out_dir: Path, base_url: str, chrom: str) -> Path:
    normalized_url = base_url.rstrip("/")
    source_id = hashlib.sha256(normalized_url.encode()).hexdigest()[:16]
    return out_dir / "_intermediate" / "clusters" / source_id / chrom


def prepare_cluster_files(
    chrom: str,
    base_url: str,
    cache_dir: Path,
    cluster_dir: Path | None,
    workers: int,
    retries: int,
) -> tuple[list[Path], dict]:
    if cluster_dir is not None:
        paths = sorted(cluster_dir.glob(f"{chrom}_*-*_cluster.tsv"), key=parse_cluster_filename)
        if not paths:
            raise ValueError(f"No {chrom} cluster files found in {cluster_dir}")
        return paths, {
            "source": "local",
            "directory": str(cluster_dir.resolve()),
            "discovered": len(paths),
            "downloaded": 0,
            "reused": len(paths),
        }

    names = discover_cluster_filenames(base_url, chrom, retries)
    cache_dir.mkdir(parents=True, exist_ok=True)
    reused = [cache_dir / name for name in names if (cache_dir / name).exists()]
    pending = [name for name in names if not (cache_dir / name).exists()]
    failures: list[str] = []
    log.info("Preparing %d cluster files (%d cached, %d to download)", len(names), len(reused), len(pending))

    def download(name: str) -> Path:
        url = f"{base_url.rstrip('/')}/{chrom}/clusters/{name}"
        destination = cache_dir / name
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        temporary.write_bytes(request_with_retries(url, retries))
        temporary.replace(destination)
        return destination

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(download, name): name for name in pending}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:  # individual failures are collected and reported together
                name = futures[future]
                log.error("Failed to download %s: %s", name, exc)
                failures.append(name)
    if failures:
        raise RuntimeError(
            f"Failed to download {len(failures)} cluster file(s); rerun to resume: {sorted(failures)}"
        )
    paths = sorted((cache_dir / name for name in names), key=parse_cluster_filename)
    return paths, {
        "source": "download",
        "base_url": base_url.rstrip("/"),
        "cache_dir": str(cache_dir.resolve()),
        "discovered": len(names),
        "downloaded": len(pending),
        "reused": len(reused),
    }


def load_sv_table(path: Path) -> tuple[list[str], list[SVRecord]]:
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[: len(METADATA_COLUMNS)] != METADATA_COLUMNS:
            raise ValueError(f"Unexpected SV table metadata columns: {header[:len(METADATA_COLUMNS)]}")
        samples = header[len(METADATA_COLUMNS) :]
        records = []
        for row in reader:
            records.append(
                SVRecord(
                    metadata=row[: len(METADATA_COLUMNS)],
                    start=int(row[3]),
                    end=int(row[4]),
                    dosages=np.fromiter(
                        (
                            -1 if (dosage := genotype_dosage(genotype)) is None else dosage
                            for genotype in row[len(METADATA_COLUMNS) :]
                        ),
                        dtype=np.int8,
                    ),
                )
            )
    return samples, records


def parse_haplotype_id(value: str) -> tuple[str, str, int, int, int]:
    match = HAPLOTYPE_RE.fullmatch(value)
    if not match:
        raise ValueError(f"Unexpected haplotype identifier: {value}")
    return (
        canonical_sample_id(match.group("sample")),
        match.group("chrom"),
        int(match.group("start")),
        int(match.group("end")),
        int(match.group("haplotype")),
    )


def load_block_membership(
    path: Path,
    vcf_samples: set[str],
) -> tuple[BlockMembership, int, int]:
    chrom, start, end = parse_cluster_filename(path)
    haploblock_id = f"{chrom}_{start}_{end}"
    total_counts: Counter[str] = Counter()
    by_sample: dict[str, list[str | None]] = {}
    total_rows = 0
    used_rows = 0

    with path.open() as handle:
        for line in handle:
            cluster_id, member_id = line.rstrip("\n").split("\t")
            sample, member_chrom, member_start, member_end, haplotype = parse_haplotype_id(member_id)
            if (member_chrom, member_start, member_end) != (chrom, start, end):
                raise ValueError(f"Membership coordinates do not match filename in {path}: {member_id}")
            total_rows += 1
            total_counts[cluster_id] += 1
            if sample not in vcf_samples:
                continue
            used_rows += 1
            clusters = by_sample.setdefault(sample, [None, None])
            clusters[haplotype] = cluster_id

    complete_samples = {
        sample: (clusters[0], clusters[1])
        for sample, clusters in by_sample.items()
        if clusters[0] is not None and clusters[1] is not None
    }
    return (
        BlockMembership(
            haploblock_id=haploblock_id,
            chrom=chrom,
            start=start,
            end=end,
            sample_clusters=complete_samples,
            total_cluster_haplotypes=total_counts,
        ),
        total_rows,
        used_rows,
    )


@cache
def genotype_dosage(genotype: str) -> int | None:
    if "." in genotype:
        return None
    alleles = genotype.replace("|", "/").split("/")
    return sum(int(allele) != 0 for allele in alleles)


def fit_cluster_probabilities(
    cluster0: np.ndarray,
    cluster1: np.ndarray,
    dosage: np.ndarray,
    n_clusters: int,
    max_iterations: int,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray, bool, int]:
    """Fit per-cluster alternate-haplotype probabilities and return expected counts."""
    callable_counts = np.bincount(
        np.concatenate([cluster0, cluster1]), minlength=n_clusters
    ).astype(float)
    initial_values = np.concatenate([dosage / 2.0, dosage / 2.0])
    initial_clusters = np.concatenate([cluster0, cluster1])
    expected_counts = np.bincount(initial_clusters, weights=initial_values, minlength=n_clusters)
    probabilities = np.divide(
        expected_counts,
        callable_counts,
        out=np.full(n_clusters, np.nan),
        where=callable_counts > 0,
    )
    probabilities[callable_counts > 0] = np.clip(probabilities[callable_counts > 0], 1e-6, 1 - 1e-6)

    converged = False
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        expected0 = np.zeros(len(dosage), dtype=float)
        expected1 = np.zeros(len(dosage), dtype=float)
        expected0[dosage == 2] = 1.0
        expected1[dosage == 2] = 1.0
        heterozygous = dosage == 1
        p0 = probabilities[cluster0[heterozygous]]
        p1 = probabilities[cluster1[heterozygous]]
        numerator = p0 * (1 - p1)
        denominator = numerator + (1 - p0) * p1
        posterior0 = np.divide(numerator, denominator, out=np.full_like(numerator, 0.5), where=denominator > 0)
        expected0[heterozygous] = posterior0
        expected1[heterozygous] = 1 - posterior0

        expected_counts = np.bincount(cluster0, weights=expected0, minlength=n_clusters)
        expected_counts += np.bincount(cluster1, weights=expected1, minlength=n_clusters)
        updated = np.divide(
            expected_counts,
            callable_counts,
            out=np.full(n_clusters, np.nan),
            where=callable_counts > 0,
        )
        updated[callable_counts > 0] = np.clip(updated[callable_counts > 0], 1e-6, 1 - 1e-6)
        difference = np.nanmax(np.abs(updated - probabilities))
        probabilities = updated
        if difference < tolerance:
            converged = True
            break
    return probabilities, expected_counts, converged, iterations


def wilson_interval(probability: float, count: int) -> tuple[float, float]:
    if count == 0 or math.isnan(probability):
        return math.nan, math.nan
    z = 1.96
    denominator = 1 + z * z / count
    center = (probability + z * z / (2 * count)) / denominator
    margin = z * math.sqrt(probability * (1 - probability) / count + z * z / (4 * count * count)) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


FRACTION_BIN_NAMES = ["[0,0.25)", "[0.25,0.5)", "[0.5,0.75)", "[0.75,1]", "missing"]


def fraction_bin_counts(values: np.ndarray) -> np.ndarray:
    missing = np.isnan(values)
    present = values[~missing]
    return np.asarray(
        [
            np.count_nonzero(present < 0.25),
            np.count_nonzero((present >= 0.25) & (present < 0.5)),
            np.count_nonzero((present >= 0.5) & (present < 0.75)),
            np.count_nonzero(present >= 0.75),
            np.count_nonzero(missing),
        ],
        dtype=np.int64,
    )


def prepare_block(block: BlockMembership, sample_index: dict[str, int]) -> PreparedBlock:
    paired_samples = [sample for sample in block.sample_clusters if sample in sample_index]
    cluster_ids = sorted(
        {cluster for sample in paired_samples for cluster in block.sample_clusters[sample]}
    )
    cluster_index = {cluster: index for index, cluster in enumerate(cluster_ids)}
    cluster0 = np.fromiter(
        (cluster_index[block.sample_clusters[sample][0]] for sample in paired_samples),
        dtype=np.int32,
    )
    cluster1 = np.fromiter(
        (cluster_index[block.sample_clusters[sample][1]] for sample in paired_samples),
        dtype=np.int32,
    )
    represented_counts = np.bincount(
        np.concatenate([cluster0, cluster1]), minlength=len(cluster_ids)
    )
    return PreparedBlock(
        membership=block,
        sample_indices=np.fromiter(
            (sample_index[sample] for sample in paired_samples), dtype=np.int32
        ),
        cluster_ids=cluster_ids,
        cluster0=cluster0,
        cluster1=cluster1,
        represented_counts=represented_counts,
    )


def infer_sv_block(
    record: SVRecord,
    block: PreparedBlock,
    association_threshold: float,
    posterior_threshold: float,
    max_iterations: int,
    tolerance: float,
) -> dict:
    dosage_values = record.dosages[block.sample_indices]
    callable_mask = dosage_values >= 0
    dosage_array = dosage_values[callable_mask]
    cluster0_array = block.cluster0[callable_mask]
    cluster1_array = block.cluster1[callable_mask]
    n_clusters = len(block.cluster_ids)

    if len(dosage_array):
        probabilities, expected_counts, converged, iterations = fit_cluster_probabilities(
            cluster0_array,
            cluster1_array,
            dosage_array,
            n_clusters,
            max_iterations,
            tolerance,
        )
        callable_counts = np.bincount(
            np.concatenate([cluster0_array, cluster1_array]), minlength=n_clusters
        )
    else:
        probabilities = np.full(n_clusters, np.nan)
        expected_counts = np.zeros(n_clusters)
        callable_counts = np.zeros(n_clusters, dtype=int)
        converged = False
        iterations = 0

    required_callable = np.minimum(3, block.represented_counts)
    associated_mask = (
        (block.represented_counts > 0)
        & (callable_counts >= required_callable)
        & ~np.isnan(probabilities)
        & (probabilities >= association_threshold)
    )
    associated_indices = np.flatnonzero(associated_mask)
    associated_ids = [block.cluster_ids[index] for index in associated_indices]
    call_rates = np.divide(
        callable_counts,
        block.represented_counts,
        out=np.full(n_clusters, np.nan),
        where=block.represented_counts > 0,
    )
    evidence = []
    for index in associated_indices:
        cluster_id = block.cluster_ids[index]
        represented = int(block.represented_counts[index])
        callable_count = int(callable_counts[index])
        probability = float(probabilities[index])
        ci_low, ci_high = wilson_interval(probability, callable_count)
        evidence.append(
            {
                "cluster_id": cluster_id,
                "cluster_haplotypes_total": block.membership.total_cluster_haplotypes[cluster_id],
                "cluster_haplotypes_in_vcf": represented,
                "callable_haplotypes": callable_count,
                "required_callable_haplotypes": int(required_callable[index]),
                "call_rate": float(call_rates[index]),
                "expected_alt_haplotypes": float(expected_counts[index]),
                "sv_probability": probability,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "evidence_tier": "low" if callable_count < 3 else "standard",
            }
        )

    assignment_counts = np.zeros(3, dtype=np.int64)
    assignment_confidence_bins = np.zeros(len(FRACTION_BIN_NAMES), dtype=np.int64)
    heterozygous = dosage_array == 1
    heterozygous_cluster0 = cluster0_array[heterozygous]
    heterozygous_cluster1 = cluster1_array[heterozygous]
    same_cluster = heterozygous_cluster0 == heterozygous_cluster1
    assignment_counts[0] = np.count_nonzero(same_cluster)
    assignment_confidence_bins[3] = assignment_counts[0]
    different_cluster0 = heterozygous_cluster0[~same_cluster]
    different_cluster1 = heterozygous_cluster1[~same_cluster]
    if len(different_cluster0):
        p0 = probabilities[different_cluster0]
        p1 = probabilities[different_cluster1]
        numerator = p0 * (1 - p1)
        denominator = numerator + (1 - p0) * p1
        posterior0 = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, 0.5),
            where=denominator > 0,
        )
        confidence = np.maximum(posterior0, 1 - posterior0)
        assigned = confidence >= posterior_threshold
        assignment_counts[1] = np.count_nonzero(assigned)
        assignment_counts[2] = len(assigned) - assignment_counts[1]
        assignment_confidence_bins += fraction_bin_counts(confidence)

    association_class = "multiple" if len(associated_ids) > 1 else "single" if associated_ids else "zero"
    return {
        "evidence": evidence,
        "clusters_evaluated": n_clusters,
        "probability_bins": fraction_bin_counts(probabilities),
        "call_rate_bins": fraction_bin_counts(call_rates),
        "assignment_counts": assignment_counts,
        "assignment_confidence_bins": assignment_confidence_bins,
        "associated_ids": associated_ids,
        "association_class": association_class,
        "callable_samples": len(dosage_array),
        "heterozygous_samples": int(np.sum(dosage_array == 1)),
        "converged": converged,
        "iterations": iterations,
    }


def evaluate_block(task: tuple) -> list[tuple[int, dict]]:
    (
        block,
        indexed_records,
        association_threshold,
        posterior_threshold,
        max_iterations,
        tolerance,
    ) = task
    return [
        (
            record_index,
            infer_sv_block(
                record,
                block,
                association_threshold,
                posterior_threshold,
                max_iterations,
                tolerance,
            ),
        )
        for record_index, record in indexed_records
    ]


def write_results(
    chrom: str,
    sv_table: Path,
    cluster_paths: list[Path],
    out_dir: Path,
    association_threshold: float,
    posterior_threshold: float,
    max_iterations: int,
    tolerance: float,
    download_qc: dict,
    workers: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_dir = out_dir / "debug_and_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    samples, records = load_sv_table(sv_table)
    sample_set = set(samples)
    log.info("Loaded %d SVs across %d samples; reading cluster memberships", len(records), len(samples))

    blocks = []
    membership_rows = 0
    membership_rows_used = 0
    for path in cluster_paths:
        block, total_rows, used_rows = load_block_membership(path, sample_set)
        blocks.append(block)
        membership_rows += total_rows
        membership_rows_used += used_rows

    log.info("Loaded %d cluster blocks; evaluating overlapping SVs", len(blocks))
    sample_index = {sample: index for index, sample in enumerate(samples)}
    prepared_blocks = [prepare_block(block, sample_index) for block in blocks]

    haploblock_path = out_dir / f"haploblocks.{chrom}.tsv"
    with haploblock_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["haploblock_id", "chrom", "start", "end"])
        for block in blocks:
            writer.writerow([block.haploblock_id, block.chrom, block.start, block.end])

    cluster_membership_path = out_dir / f"cluster_memberships.{chrom}.tsv"
    with cluster_membership_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["haploblock_id", "chrom", "start", "end", "sample_id", "haplotype", "cluster_id"])
        for block in blocks:
            for sample_id, cluster_ids in sorted(block.sample_clusters.items()):
                for haplotype, cluster_id in enumerate(cluster_ids):
                    writer.writerow(
                        [block.haploblock_id, block.chrom, block.start, block.end, sample_id, haplotype, cluster_id]
                    )

    starts = np.asarray([record.start for record in records])
    ends = np.asarray([record.end for record in records])
    overlap_counts = np.zeros(len(records), dtype=int)
    block_record_indices = []
    for block in blocks:
        indices = np.flatnonzero((starts < block.end) & (ends > block.start))
        block_record_indices.append(indices)
        overlap_counts[indices] += 1

    downstream_fields = [
        *METADATA_COLUMNS, "haploblock_id", "block_start", "block_end", "overlaps_multiple_blocks",
        "cluster_id", "cluster_haplotypes_total", "cluster_haplotypes_in_vcf", "callable_haplotypes",
        "required_callable_haplotypes", "call_rate", "expected_alt_haplotypes", "sv_probability",
        "ci95_low", "ci95_high", "evidence_tier", "model_converged", "model_iterations",
    ]
    summary_fields = [
        *METADATA_COLUMNS, "haploblock_id", "block_start", "block_end", "overlaps_multiple_blocks",
        "callable_samples", "heterozygous_samples", "clusters_evaluated", "associated_clusters",
        "associated_cluster_ids", "association_class", "model_converged", "model_iterations",
    ]
    association_counts: Counter[str] = Counter()
    converged_models = 0
    nonconverged_models = 0
    no_callable_pairs = 0
    sv_block_pairs = 0
    cluster_evaluations = 0
    associated_cluster_rows = 0
    low_evidence_associations = 0
    probability_bin_counts = np.zeros(len(FRACTION_BIN_NAMES), dtype=np.int64)
    call_rate_bin_counts = np.zeros(len(FRACTION_BIN_NAMES), dtype=np.int64)
    assignment_status_counts = np.zeros(3, dtype=np.int64)
    assignment_confidence_bin_counts = np.zeros(len(FRACTION_BIN_NAMES), dtype=np.int64)

    block_tasks = [
        (
            block,
            [(int(index), records[int(index)]) for index in record_indices],
            association_threshold,
            posterior_threshold,
            max_iterations,
            tolerance,
        )
        for block, record_indices in zip(prepared_blocks, block_record_indices)
    ]

    with (
        (out_dir / f"sv_to_clusters.{chrom}.tsv").open("w", newline="") as downstream_handle,
        (out_dir / f"sv_block_summary.{chrom}.tsv").open("w", newline="") as summary_handle,
    ):
        downstream_writer = csv.DictWriter(downstream_handle, fieldnames=downstream_fields, delimiter="\t", lineterminator="\n")
        summary_writer = csv.DictWriter(summary_handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        downstream_writer.writeheader()
        summary_writer.writeheader()

        if workers > 1 and len(block_tasks) > 1:
            try:
                executor = ProcessPoolExecutor(
                    max_workers=min(workers, len(block_tasks)),
                    mp_context=multiprocessing.get_context("spawn"),
                )
                analysis_workers = min(workers, len(block_tasks))
                evaluated_blocks = executor.map(
                    evaluate_block,
                    block_tasks,
                    chunksize=max(1, len(block_tasks) // (workers * 4)),
                )
            except (NotImplementedError, PermissionError) as error:
                log.warning("Process workers unavailable (%s); evaluating blocks sequentially", error)
                executor = None
                analysis_workers = 1
                evaluated_blocks = map(evaluate_block, block_tasks)
        else:
            executor = None
            analysis_workers = 1
            evaluated_blocks = map(evaluate_block, block_tasks)

        try:
            for prepared_block, evaluated_pairs in zip(prepared_blocks, evaluated_blocks):
                block = prepared_block.membership
                for record_index, result in evaluated_pairs:
                    record = records[record_index]
                    sv_block_pairs += 1
                    association_counts[result["association_class"]] += 1
                    if result["callable_samples"] == 0:
                        no_callable_pairs += 1
                    elif result["converged"]:
                        converged_models += 1
                    else:
                        nonconverged_models += 1
                    common = dict(zip(METADATA_COLUMNS, record.metadata))
                    common.update(
                        {
                            "haploblock_id": block.haploblock_id,
                            "block_start": block.start,
                            "block_end": block.end,
                            "overlaps_multiple_blocks": overlap_counts[record_index] > 1,
                        }
                    )
                    cluster_evaluations += result["clusters_evaluated"]
                    probability_bin_counts += result["probability_bins"]
                    call_rate_bin_counts += result["call_rate_bins"]
                    for evidence in result["evidence"]:
                        associated_cluster_rows += 1
                        low_evidence_associations += evidence["evidence_tier"] == "low"
                        downstream_writer.writerow(
                            {
                                **common,
                                **{field: evidence[field] for field in downstream_fields if field in evidence},
                                "model_converged": result["converged"],
                                "model_iterations": result["iterations"],
                            }
                        )
                    assignment_status_counts += result["assignment_counts"]
                    assignment_confidence_bin_counts += result["assignment_confidence_bins"]
                    summary_writer.writerow(
                        {
                            **common,
                            "callable_samples": result["callable_samples"],
                            "heterozygous_samples": result["heterozygous_samples"],
                            "clusters_evaluated": result["clusters_evaluated"],
                            "associated_clusters": len(result["associated_ids"]),
                            "associated_cluster_ids": ",".join(result["associated_ids"]),
                            "association_class": result["association_class"],
                            "model_converged": result["converged"],
                            "model_iterations": result["iterations"],
                        }
                    )
        finally:
            if executor is not None:
                executor.shutdown()

    analyzable_pairs = sv_block_pairs - no_callable_pairs
    qc = {
        "chrom": chrom,
        "downstream_output": str((out_dir / f"sv_to_clusters.{chrom}.tsv").resolve()),
        "haploblock_output": str(haploblock_path.resolve()),
        "cluster_membership_output": str(cluster_membership_path.resolve()),
        "sv_block_summary_output": str((out_dir / f"sv_block_summary.{chrom}.tsv").resolve()),
        "sv_table": str(sv_table.resolve()),
        "sv_records": len(records),
        "vcf_samples": len(samples),
        "cluster_download": download_qc,
        "cluster_files": len(blocks),
        "cluster_membership_rows": membership_rows,
        "cluster_membership_rows_matching_vcf": membership_rows_used,
        "sv_block_pairs": sv_block_pairs,
        "sv_block_pairs_without_callable_genotypes": no_callable_pairs,
        "converged_models": converged_models,
        "nonconverged_models": nonconverged_models,
        "cluster_evaluations": cluster_evaluations,
        "associated_cluster_rows": associated_cluster_rows,
        "low_evidence_associations": low_evidence_associations,
        "cluster_probability_bins": {
            name: int(count)
            for name, count in zip(FRACTION_BIN_NAMES, probability_bin_counts)
            if count
        },
        "cluster_call_rate_bins": {
            name: int(count)
            for name, count in zip(FRACTION_BIN_NAMES, call_rate_bin_counts)
            if count
        },
        "heterozygote_assignment_status_counts": {
            name: int(count)
            for name, count in zip(["same_cluster", "assigned", "ambiguous"], assignment_status_counts)
            if count
        },
        "heterozygote_assignment_confidence_bins": {
            name: int(count)
            for name, count in zip(FRACTION_BIN_NAMES, assignment_confidence_bin_counts)
            if count
        },
        "association_class_counts": dict(association_counts),
        "multiple_cluster_fraction": association_counts["multiple"] / analyzable_pairs if analyzable_pairs else None,
        "association_threshold": association_threshold,
        "posterior_assignment_threshold": posterior_threshold,
        "max_iterations": max_iterations,
        "tolerance": tolerance,
        "cluster_cache_retained": download_qc["source"] == "download",
        "analysis_workers": analysis_workers,
    }
    (qc_dir / f"method_qc.{chrom}.json").write_text(json.dumps(qc, indent=2) + "\n")
    log.info("Evaluated %d SV-block pairs across %d cluster files", sv_block_pairs, len(blocks))


def process_chromosome(
    chrom: str,
    sv_table: Path,
    out_dir: Path,
    cluster_base_url: str,
    cluster_root: Path | None,
    workers: int,
    retries: int,
    association_threshold: float,
    posterior_threshold: float,
    max_iterations: int,
    tolerance: float,
) -> None:
    cluster_dir = None
    if cluster_root is not None:
        per_chrom_dir = cluster_root / chrom
        cluster_dir = per_chrom_dir if per_chrom_dir.exists() else cluster_root
    cluster_cache_dir = cluster_cache_directory(out_dir, cluster_base_url, chrom)
    cluster_paths, download_qc = prepare_cluster_files(
        chrom,
        cluster_base_url,
        cluster_cache_dir,
        cluster_dir,
        workers,
        retries,
    )
    write_results(
        chrom,
        sv_table,
        cluster_paths,
        out_dir,
        association_threshold,
        posterior_threshold,
        max_iterations,
        tolerance,
        download_qc,
        workers,
    )



def remove_legacy_outputs(out_dir: Path) -> None:
    for path in (
        out_dir / "sample_id_map.tsv",
        out_dir / "cluster_evidence.tsv",
        out_dir / "sv_cluster_associations.tsv",
        out_dir / "heterozygote_assignments.tsv",
        out_dir / "sv_block_summary.tsv",
        out_dir / "cluster_qc.json",
        out_dir / "sv_to_clusters.tsv",
        out_dir / "debug_and_qc" / "sv_block_qc.tsv",
        out_dir / "debug_and_qc" / "method_qc.json",
        out_dir / "debug_and_qc" / "vcf_qc.json",
    ):
        path.unlink(missing_ok=True)


def chromosome_is_complete(
    chrom: str,
    vcf_path: Path,
    out_dir: Path,
    cluster_base_url: str,
    cluster_root: Path | None,
    max_sv_id_length: int,
    association_threshold: float,
    posterior_threshold: float,
    max_iterations: int,
    tolerance: float,
) -> bool:
    qc_dir = out_dir / "debug_and_qc"
    required = [
        out_dir / "samples.tsv",
        out_dir / f"sv_genotypes.{chrom}.tsv",
        out_dir / f"sv_to_clusters.{chrom}.tsv",
        out_dir / f"haploblocks.{chrom}.tsv",
        out_dir / f"cluster_memberships.{chrom}.tsv",
        out_dir / f"sv_block_summary.{chrom}.tsv",
        qc_dir / f"method_qc.{chrom}.json",
        qc_dir / f"vcf_qc.{chrom}.json",
    ]
    if not all(path.exists() for path in required):
        return False
    vcf_qc = json.loads((qc_dir / f"vcf_qc.{chrom}.json").read_text())
    method_qc = json.loads((qc_dir / f"method_qc.{chrom}.json").read_text())
    cluster_download = method_qc.get("cluster_download", {})
    if cluster_root is None:
        cluster_source_matches = (
            cluster_download.get("source") == "download"
            and cluster_download.get("base_url") == cluster_base_url.rstrip("/")
        )
    else:
        per_chrom_dir = cluster_root / chrom
        cluster_dir = per_chrom_dir if per_chrom_dir.exists() else cluster_root
        cluster_source_matches = (
            cluster_download.get("source") == "local"
            and cluster_download.get("directory") == str(cluster_dir.resolve())
        )
    return (
        vcf_qc.get("vcf") == str(vcf_path.resolve())
        and vcf_qc.get("sv_id_max_length") == max_sv_id_length
        and cluster_source_matches
        and method_qc.get("association_threshold") == association_threshold
        and method_qc.get("posterior_assignment_threshold") == posterior_threshold
        and method_qc.get("max_iterations") == max_iterations
        and method_qc.get("tolerance") == tolerance
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, default=Path("input/1kgp_ont_cohort.postfilter.full.vcf.gz"))
    parser.add_argument("--chroms", default="all", help="Comma-separated chromosomes or 'all' for chr1-22,X")
    parser.add_argument("--cluster-base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--cluster-root", type=Path, default=None, help="Use local <root>/<chrom> cluster files")
    parser.add_argument("--out-dir", type=Path, default=Path("stage1_output"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--max-sv-id-length", type=int, default=80)
    parser.add_argument("--association-threshold", type=float, default=0.75)
    parser.add_argument("--posterior-threshold", type=float, default=0.75)
    parser.add_argument("--max-iterations", type=int, default=25)
    parser.add_argument("--tolerance", type=float, default=1e-5)
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> None:
    chroms = (
        ALL_CHROMS
        if args.chroms.strip().lower() == "all"
        else [normalize_chrom(chrom.strip()) for chrom in args.chroms.split(",") if chrom.strip()]
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pending_chroms = [
        chrom
        for chrom in chroms
        if not chromosome_is_complete(
            chrom,
            args.vcf,
            args.out_dir,
            args.cluster_base_url,
            args.cluster_root,
            args.max_sv_id_length,
            args.association_threshold,
            args.posterior_threshold,
            args.max_iterations,
            args.tolerance,
        )
    ]
    completed_at_start = [chrom for chrom in chroms if chrom not in pending_chroms]
    if completed_at_start:
        log.info("Reusing completed outputs for: %s", ", ".join(completed_at_start))

    chrom_workers = min(args.threads, len(pending_chroms)) if pending_chroms else 1
    workers_per_chromosome = max(1, args.threads // chrom_workers)
    sv_tables = {}
    if pending_chroms:
        sv_tables = reusable_sv_tables(
            args.vcf, pending_chroms, args.out_dir, args.max_sv_id_length
        )
        if sv_tables is None:
            sv_tables = split_vcf_by_chromosome(
                args.vcf, pending_chroms, args.out_dir, args.max_sv_id_length, chrom_workers
            )

    if pending_chroms:
        chromosome_tasks = [
            (
                chrom,
                sv_tables[chrom],
                args.out_dir,
                args.cluster_base_url,
                args.cluster_root,
                workers_per_chromosome,
                args.retries,
                args.association_threshold,
                args.posterior_threshold,
                args.max_iterations,
                args.tolerance,
            )
            for chrom in pending_chroms
        ]
        if chrom_workers == 1:
            for task in chromosome_tasks:
                process_chromosome(*task)
                log.info("Completed %s", task[0])
        else:
            with ThreadPoolExecutor(max_workers=chrom_workers) as executor:
                futures = {
                    executor.submit(process_chromosome, *task): task[0]
                    for task in chromosome_tasks
                }
                for future in as_completed(futures):
                    chrom = futures[future]
                    future.result()
                    log.info("Completed %s", chrom)

    qc_dir = args.out_dir / "debug_and_qc"
    completed_chroms = [
        chrom for chrom in ALL_CHROMS if (qc_dir / f"method_qc.{chrom}.json").exists()
    ]
    run_qc = {
        "vcf": str(args.vcf.resolve()),
        "completed_chromosomes": completed_chroms,
        "last_invocation": {
            "requested_chromosomes": chroms,
            "threads": args.threads,
            "chromosome_workers": chrom_workers,
            "download_workers_per_chromosome": workers_per_chromosome,
            "analysis_workers_per_chromosome": workers_per_chromosome,
        },
        "outputs": {
            chrom: {
                "sv_genotypes": str((args.out_dir / f"sv_genotypes.{chrom}.tsv").resolve()),
                "sv_to_clusters": str((args.out_dir / f"sv_to_clusters.{chrom}.tsv").resolve()),
                "haploblocks": str((args.out_dir / f"haploblocks.{chrom}.tsv").resolve()),
                "cluster_memberships": str((args.out_dir / f"cluster_memberships.{chrom}.tsv").resolve()),
                "sv_block_summary": str((args.out_dir / f"sv_block_summary.{chrom}.tsv").resolve()),
            }
            for chrom in completed_chroms
        },
    }
    (qc_dir / "run_qc.json").write_text(json.dumps(run_qc, indent=2) + "\n")
    if completed_chroms == ALL_CHROMS:
        remove_legacy_outputs(args.out_dir)


def main(argv: list[str] | None = None) -> None:
    run(parse_args(argv))


if __name__ == "__main__":
    main()
