"""Stage 0: Data ingestion & harmonization for the Haploblock_SV pipeline.

Loads (or synthesizes, as a fallback) three inputs:
  1. dbVar nstd152 structural variant calls (VCF)
  2. haploblock coordinates/metadata + population-cluster labels (data.haploblocks.org)
  3. 1000 Genomes sample -> superpopulation metadata

...checks they share a genome build (lifting over with pyliftover if a chain
file is supplied and they don't), then writes everything plus a shared
config.yaml that every later pipeline stage reads. Any source that is not
provided, unreachable, or fails to parse falls back to a synthetic
stand-in so downstream stages are never blocked.

Coordinate convention: all `start`/`end` columns in the output tables are
0-based, half-open (BED-style), matching the synthetic generator and a
BED-style haploblock source. VCF POS (1-based) is converted to 0-based on
load. Genome-build harmonization only applies to the SV and haploblock
tables, which carry genomic coordinates -- the sample->superpopulation
panel has none, so it is never lifted.

--------------------------------------------------------------------------
Contract for pipeline/stage1_qc.py (its input; see stage1_qc.py itself
for what it produces)
--------------------------------------------------------------------------
Stage 1 locates its inputs via config.yaml -- read
config["paths"]["sv_calls"], config["paths"]["haploblocks"], and
config["paths"]["sample_metadata"] (each an absolute path) rather than
hardcoding "example_data/..." paths, since --out-dir can move them.

All three table files are TSV (tab-separated, header row via
pandas.to_csv(sep="\t")). Coordinates are 0-based, half-open; chrom is
always UCSC-style ("chr1".."chr22", "chrX").

  sv_calls.tsv        sv_id (str), chrom (str), start (int), end (int),
                       sv_type (str, one of DEL/DUP/INV/INS), imprecise
                       (bool -- dbVar never sets VCF FILTER, so this is
                       nstd152's closest thing to a confidence flag),
                       length (float, nullable -- end-start for DEL/DUP/INV;
                       for INS this is NOT end-start (always ~1bp, an
                       insertion's reference *position*, not its size) but
                       a length recovered from supporting calls' SVLEN,
                       NaN if no call resolved one), then one column per
                       sample -- column names match sample_metadata.tsv's
                       sample_id values exactly (join key, not a fixed set
                       of names). The per-sample cell is a 0/1/2 dosage int
                       for the dbVar and synthetic paths, but the raw
                       FORMAT/GT string ("0|1", "1/1", "./.") when the
                       calls came from --vcf (a standard multi-sample VCF);
                       downstream stages must handle both (a GT string that
                       is all-numeric-and-separators vs. a bare int).

  haploblocks.tsv     haploblock_id (str), chrom (str), start (int),
                       end (int), n_snps (float, always NaN -- not published
                       by data.haploblocks.org's boundary/hash files), n_clusters
                       (float, always NaN -- would need the thousands of
                       per-block chrN_cluster_hashes_<start>-<end>.tsv files,
                       not fetched by default; see fetch_haploblock_boundaries()),
                       hash_length (float, NaN if a chromosome's
                       *_haploblock_hashes.tsv wasn't available/parseable),
                       cluster_diff_score (float, always NaN for real
                       data -- an actual differentiation score needs
                       population labels and belongs in Stage 6, not here).

  sample_metadata.tsv sample_id (str), population (str -- the fine-grained
                       label downstream stages group by; a 1000G population
                       code like CHS/PUR/YRI where known, else UNKNOWN),
                       superpopulation (str, one of AFR/AMR/EAS/EUR/SAS, or
                       UNKNOWN). Stages 4/6/7 read `population`, never a
                       hardcoded list of the five superpopulations.

config.yaml also carries values Stage 1 (and later stages) should read
rather than re-hardcode: genome_build (str), thresholds.af_common_threshold,
thresholds.boundary_distance_bp, thresholds.min_sv_count_per_block,
thresholds.min_sv_length / max_sv_length (bp, Stage 1's size filter),
thresholds.drop_imprecise (bool, Stage 1's confidence filter), and
seeds.permutation_seed / seeds.umap_seed.
"""

from __future__ import annotations

import argparse
import logging
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage0_ingest")

SUPERPOPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SV_TYPES = ["DEL", "DUP", "INV", "INS"]
SV_TYPE_WEIGHTS = [0.35, 0.25, 0.15, 0.25]
CHROMS = ["chr1", "chr2", "chr3"]

DBVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf"
HAPLOBLOCK_HASH_BASE_URL = "https://data.haploblocks.org/haploblock_hashes/1000G"
ALL_1000G_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# nstd152's variant_call.vcf.gz gives sample IDs (INFO/SAMPLE) but no population
# labels. As of 2026-08 the study's calls cover exactly 9 samples from 3 known
# HGSVC trios; hardcoded here since there is no --panel-source lookup for them.
# If a future dbVar release adds samples outside this set, they surface as
# UNKNOWN (see the warning in main()) rather than being silently mislabeled.
HGSVC_TRIO_SUPERPOPULATIONS = {
    "HG00512": "EAS", "HG00513": "EAS", "HG00514": "EAS",  # CHS trio
    "HG00731": "AMR", "HG00732": "AMR", "HG00733": "AMR",  # PUR trio
    "NA19238": "AFR", "NA19239": "AFR", "NA19240": "AFR",  # YRI trio
}
HGSVC_TRIO_POPULATIONS = {
    "HG00512": "CHS", "HG00513": "CHS", "HG00514": "CHS",
    "HG00731": "PUR", "HG00732": "PUR", "HG00733": "PUR",
    "NA19238": "YRI", "NA19239": "YRI", "NA19240": "YRI",
}


# --------------------------------------------------------------------------
# Fetching real sources (best-effort; any failure falls back to synthetic)
# --------------------------------------------------------------------------

def fetch_or_load(source: str | None, dest_dir: Path, label: str) -> Path | None:
    """Return a local path for `source` (URL or local path), or None on failure."""
    if not source:
        return None
    if source.startswith("http://") or source.startswith("https://"):
        try:
            import requests

            dest = dest_dir / f"raw_{label}{Path(source).suffix or '.dat'}"
            log.info("Downloading %s from %s", label, source)
            resp = requests.get(source, timeout=30)
            resp.raise_for_status()
            dest.write_bytes(resp.content)
            return dest
        except Exception as exc:  # noqa: BLE001 - best-effort fetch, any failure -> fallback
            log.warning("Failed to download %s (%s): %s", label, source, exc)
            return None
    path = Path(source)
    if not path.exists():
        log.warning("Local path for %s does not exist: %s", label, path)
        return None
    return path


def download_file(url: str, dest_path: Path, timeout: int = 60) -> bool:
    """Best-effort single-file download. Returns True on success, logs and returns False otherwise."""
    try:
        import requests

        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        dest_path.write_bytes(resp.content)
        return True
    except Exception as exc:  # noqa: BLE001 - best-effort fetch, any failure -> caller falls back
        log.warning("Failed to download %s: %s", url, exc)
        return False


def fetch_dbvar_study(base_url: str, study: str, build: str, dest_dir: Path, timeout: int = 120) -> dict | None:
    """Download dbVar's paired call+region VCFs (+ .tbi indexes) for one study/build.

    Returns a dict of the four downloaded paths, or None if either VCF (the two
    files actually needed to parse SV calls) could not be fetched. A missing
    .tbi index is only a warning -- Stage 0 does a full linear parse and never
    needs random access.
    """
    names = {
        "call_vcf": f"{study}.{build}.variant_call.vcf.gz",
        "call_tbi": f"{study}.{build}.variant_call.vcf.gz.tbi",
        "region_vcf": f"{study}.{build}.variant_region.vcf.gz",
        "region_tbi": f"{study}.{build}.variant_region.vcf.gz.tbi",
    }
    paths = {}
    for label, fname in names.items():
        dest = dest_dir / fname
        if download_file(f"{base_url}/{fname}", dest, timeout=timeout):
            paths[label] = dest
            log.info("Downloaded %s (%d bytes)", fname, dest.stat().st_size)
        elif label.endswith("_tbi"):
            log.warning("Could not download index %s -- not required for a full linear parse, continuing", fname)
        else:
            log.warning("Could not download required dbVar file %s; aborting dbVar fetch", fname)
            return None
    return paths


def load_dbvar_sv_calls(call_vcf_path: Path, region_vcf_path: Path) -> pd.DataFrame | None:
    """Parse dbVar's paired call+region VCFs into one SV table with per-sample dosage columns.

    dbVar's variant_call.vcf.gz has one row per (sample, supporting call) via
    INFO/SAMPLE and INFO/REGIONID -- there is no FORMAT/GT genotype block, unlike
    a typical multi-sample VCF. variant_region.vcf.gz holds one merged, canonical
    interval + SVTYPE per REGIONID accession; we treat it as authoritative for
    coordinates/type, and use the call file only to derive per-sample dosage
    (capped at 2, i.e. treating >=2 supporting calls for the same sample+region
    as homozygous).

    Also emits `imprecise` (bool, from INFO/IMPRECISE) and `length` (float,
    nullable). dbVar never sets VCF QUAL/FILTER (both are always "."), so
    `imprecise` is the closest thing this source has to a call-confidence flag.
    `length` is `end - start` for DEL/DUP/INV, where that's meaningful, but
    variant_region.vcf.gz's END/SVLEN for INS describe an imprecise call's
    *locus* span, not the inserted sequence length, and INS's own SVLEN is
    always "." at the region level -- every usable INS length in this dataset
    only exists per-call (INFO/SVLEN in variant_call.vcf.gz), so it's collected
    there instead and folded in as each region's median across its calls.
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        log.warning("cyvcf2 not installed; cannot parse dbVar VCFs")
        return None
    try:
        region_rows = []
        n_region_skipped = 0
        for variant in VCF(str(region_vcf_path)):
            svtype = variant.INFO.get("SVTYPE")
            if svtype not in SV_TYPES:
                n_region_skipped += 1
                continue
            start = variant.POS - 1  # VCF POS is 1-based; convert to 0-based half-open
            end = variant.INFO.get("END", variant.POS)
            region_rows.append(
                {
                    "sv_id": variant.ID,
                    "chrom": variant.CHROM,
                    "start": start,
                    "end": end,
                    "sv_type": svtype,
                    "imprecise": bool(variant.INFO.get("IMPRECISE", False)),
                    "length": np.nan if svtype == "INS" else abs(end - start),
                }
            )
        if n_region_skipped:
            log.warning(
                "Skipped %d variant_region record(s) with SVTYPE outside %s (e.g. CNV)",
                n_region_skipped, SV_TYPES,
            )
        if not region_rows:
            log.warning("No usable variant_region records parsed from %s", region_vcf_path)
            return None
        region_df = pd.DataFrame(region_rows).drop_duplicates(subset="sv_id").set_index("sv_id")

        dosage_counts: dict[tuple[str, str], int] = {}
        ins_lengths: dict[str, list[int]] = {}
        n_calls_used, n_calls_skipped = 0, 0
        for variant in VCF(str(call_vcf_path)):
            sample = variant.INFO.get("SAMPLE")
            regionid_field = variant.INFO.get("REGIONID")
            if not sample or not regionid_field:
                n_calls_skipped += 1  # e.g. SAMPLESET-based cohort calls, no single sample id
                continue
            svlen_field = variant.INFO.get("SVLEN")
            for regionid in str(regionid_field).split(","):
                if regionid not in region_df.index:
                    continue  # region filtered out above (e.g. CNV) or not in this region file
                key = (regionid, sample)
                dosage_counts[key] = min(2, dosage_counts.get(key, 0) + 1)
                n_calls_used += 1
                if region_df.at[regionid, "sv_type"] == "INS" and svlen_field not in (None, "."):
                    try:
                        ins_lengths.setdefault(regionid, []).append(abs(int(svlen_field)))
                    except (TypeError, ValueError):
                        pass
        if n_calls_skipped:
            log.warning("Skipped %d variant_call record(s) with no single SAMPLE/REGIONID", n_calls_skipped)

        if ins_lengths:
            median_ins_length = {rid: float(np.median(lens)) for rid, lens in ins_lengths.items()}
            region_df["length"] = region_df["length"].fillna(pd.Series(median_ins_length))
        n_ins = (region_df["sv_type"] == "INS").sum()
        n_ins_unresolved = int(((region_df["sv_type"] == "INS") & region_df["length"].isna()).sum())
        if n_ins_unresolved:
            log.warning(
                "%d/%d INS region(s) have no resolvable length from any supporting call's SVLEN "
                "-- left as NaN, exempted from Stage 1's size filter rather than dropped",
                n_ins_unresolved, n_ins,
            )

        samples = sorted({s for _, s in dosage_counts})
        log.info("dbVar call file: %d usable (sample, call) observation(s) across %d sample(s)", n_calls_used, len(samples))
        if len(samples) < len(SUPERPOPULATIONS):
            log.warning(
                "Only %d sample(s) observed in dbVar calls (%s) -- population-specific SV "
                "classification will be severely underpowered and cannot cover all 5 superpopulations",
                len(samples), ", ".join(samples),
            )

        if dosage_counts:
            dosage_series = pd.Series(dosage_counts)
            dosage_series.index = pd.MultiIndex.from_tuples(dosage_series.index, names=["sv_id", "sample_id"])
            dosage_matrix = dosage_series.unstack("sample_id", fill_value=0).reindex(region_df.index, fill_value=0).astype(int)
        else:
            dosage_matrix = pd.DataFrame(index=region_df.index)

        return region_df.join(dosage_matrix).reset_index()
    except Exception as exc:  # noqa: BLE001
        log.warning("Failed to parse dbVar VCFs: %s", exc)
        return None


def fetch_haploblock_boundaries(base_url: str, chroms: list[str], timeout: int = 30) -> pd.DataFrame | None:
    """Fetch haploblock boundaries (and, where available, a per-block hash
    length) from data.haploblocks.org's per-chromosome aggregate files.

    Each chromosome subdirectory publishes a single
    `<chrom>_haploblock_boundaries_<chrom>.tsv` (columns START, END -- one
    row per haploblock) and usually a single `<chrom>_haploblock_hashes.tsv`
    (columns START, END, HASH -- one representative hash per haploblock,
    not the full per-cluster enumeration the individual
    `chrN_cluster_hashes_<start>-<end>.tsv` files under the same directory
    give). Fetching those thousands of individual per-block files does not
    scale to a genome-wide default (tens of thousands of HTTP requests, ~2-3
    minutes observed against the real site for a 50-blocks/chromosome cap
    alone); these two per-chromosome aggregates give one small request each
    instead (23 total for the whole genome), so `n_clusters` (which needs
    the per-block cluster enumeration) is not derivable this way and stays
    NaN -- only `hash_length` is recovered, from the single hash's length.
    """
    import requests

    frames = []
    n_chrom_failed = 0
    for chrom in chroms:
        boundaries_url = f"{base_url}/{chrom}/{chrom}_haploblock_boundaries_{chrom}.tsv"
        try:
            resp = requests.get(boundaries_url, timeout=timeout)
            resp.raise_for_status()
            boundaries = pd.read_csv(StringIO(resp.text), sep="\t")
        except Exception as exc:  # noqa: BLE001
            log.warning("Failed to fetch haploblock boundaries for %s: %s", chrom, exc)
            n_chrom_failed += 1
            continue
        start_col = _find_column(list(boundaries.columns), ["start"])
        end_col = _find_column(list(boundaries.columns), ["end"])
        if not (start_col and end_col):
            log.warning("%s: haploblock boundaries file has no START/END columns, skipping", chrom)
            n_chrom_failed += 1
            continue
        chrom_df = pd.DataFrame(
            {"start": boundaries[start_col].astype(int), "end": boundaries[end_col].astype(int)}
        )

        hashes_url = f"{base_url}/{chrom}/{chrom}_haploblock_hashes.tsv"
        try:
            resp = requests.get(hashes_url, timeout=timeout)
            resp.raise_for_status()
            # dtype=str: the HASH column is a zero-padded digit string (e.g.
            # "00000000000000000001") that pandas would otherwise silently
            # parse as an integer, stripping the leading zeros and corrupting
            # its length
            hashes = pd.read_csv(StringIO(resp.text), sep="\t", dtype=str)
            h_start_col = _find_column(list(hashes.columns), ["start"])
            h_end_col = _find_column(list(hashes.columns), ["end"])
            h_hash_col = _find_column(list(hashes.columns), ["hash"])
            if h_start_col and h_end_col and h_hash_col:
                hashes = hashes.rename(columns={h_start_col: "start", h_end_col: "end", h_hash_col: "hash_value"})
                hashes["start"] = hashes["start"].astype(int)
                hashes["end"] = hashes["end"].astype(int)
                chrom_df = chrom_df.merge(hashes[["start", "end", "hash_value"]], on=["start", "end"], how="left")
                chrom_df["hash_length"] = chrom_df["hash_value"].dropna().astype(str).str.len()
                chrom_df = chrom_df.drop(columns=["hash_value"])
            else:
                chrom_df["hash_length"] = np.nan
        except Exception as exc:  # noqa: BLE001
            log.info("No/unusable haploblock_hashes file for %s (%s); hash_length left NaN", chrom, exc)
            chrom_df["hash_length"] = np.nan

        chrom_df["chrom"] = chrom
        chrom_df["haploblock_id"] = [f"{chrom}_{s}_{e}" for s, e in zip(chrom_df["start"], chrom_df["end"])]
        chrom_df["n_snps"] = np.nan
        chrom_df["n_clusters"] = np.nan
        chrom_df["cluster_diff_score"] = np.nan
        log.info("%s: fetched %d haploblock boundaries", chrom, len(chrom_df))
        frames.append(chrom_df[["haploblock_id", "chrom", "start", "end", "n_snps", "n_clusters", "hash_length", "cluster_diff_score"]])

    if n_chrom_failed:
        log.warning("%d/%d chromosome(s) failed to fetch haploblock boundaries", n_chrom_failed, len(chroms))
    if not frames:
        return None
    result = pd.concat(frames, ignore_index=True)
    log.info("Fetched %d haploblock boundaries across %d chromosome(s)", len(result), len(frames))
    return result


def load_real_sv_calls(vcf_path: Path) -> pd.DataFrame | None:
    try:
        from cyvcf2 import VCF
    except ImportError:
        log.warning("cyvcf2 not installed; cannot parse real SV VCF %s", vcf_path)
        return None
    try:
        vcf = VCF(str(vcf_path))
        sample_ids = list(vcf.samples)
        rows = []
        n_skipped = 0
        for variant in vcf:
            svtype = variant.INFO.get("SVTYPE")
            if svtype not in SV_TYPES:
                n_skipped += 1
                continue
            start = variant.POS - 1  # VCF POS is 1-based; convert to 0-based half-open
            end = variant.INFO.get("END", variant.POS)  # END is already the 0-based-exclusive stop
            svlen_field = variant.INFO.get("SVLEN")
            try:
                length = abs(int(svlen_field)) if svlen_field not in (None, ".") else abs(end - start)
            except (TypeError, ValueError):
                length = abs(end - start)
            row = {
                "sv_id": variant.ID or f"{variant.CHROM}_{variant.POS}_{svtype}",
                "chrom": variant.CHROM,
                "start": start,
                "end": end,
                "sv_type": svtype,
                "imprecise": bool(variant.INFO.get("IMPRECISE", False)),
                "length": length if length > 0 else np.nan,
            }
            if sample_ids and variant.gt_types is not None:
                # cyvcf2 gt_types: 0=HOM_REF,1=HET,2=HOM_ALT,3=UNKNOWN -> dosage
                dosage_map = {0: 0, 1: 1, 2: 2, 3: 0}
                for sample_id, gt in zip(sample_ids, variant.gt_types):
                    row[sample_id] = dosage_map.get(int(gt), 0)
            rows.append(row)
        if n_skipped:
            log.warning("Skipped %d VCF records with unrecognized/missing SVTYPE", n_skipped)
        if not rows:
            log.warning("No usable SV records parsed from %s", vcf_path)
            return None
        return pd.DataFrame(rows)
    except Exception as exc:  # noqa: BLE001
        log.warning("Failed to parse SV VCF %s: %s", vcf_path, exc)
        return None


def load_vcf_sv_calls(vcf_path: Path) -> pd.DataFrame | None:
    """Parse a standard multi-sample SV VCF, keeping the raw FORMAT/GT string
    ("0|1", "1/1", "./.") in one column per sample.

    This is the --vcf path. Unlike load_real_sv_calls() (which collapses each
    genotype to a 0/1/2 dosage), the literal GT is preserved so phasing and
    missingness survive into sv_calls.tsv for Stage 4's per-population allele
    frequencies. Records whose INFO/SVTYPE is not one of DEL/DUP/INV/INS are
    skipped. Returns None (caller falls back) on any parse failure.
    """
    try:
        from cyvcf2 import VCF
    except ImportError:
        log.warning("cyvcf2 not installed; cannot parse --vcf %s", vcf_path)
        return None
    try:
        vcf = VCF(str(vcf_path))
        sample_ids = list(vcf.samples)
        if not sample_ids:
            log.warning("--vcf %s declares no samples in its header", vcf_path)
            return None
        rows = []
        n_skipped = 0
        for variant in vcf:
            svtype = variant.INFO.get("SVTYPE")
            if svtype not in SV_TYPES:
                n_skipped += 1
                continue
            start = variant.POS - 1  # VCF POS is 1-based; convert to 0-based half-open
            end = variant.INFO.get("END", variant.POS)
            svlen_field = variant.INFO.get("SVLEN")
            try:
                length = abs(int(svlen_field)) if svlen_field not in (None, ".") else abs(end - start)
            except (TypeError, ValueError):
                length = abs(end - start)
            row = {
                "sv_id": variant.ID or f"{variant.CHROM}_{variant.POS}_{svtype}",
                "chrom": variant.CHROM,
                "start": start,
                "end": end,
                "sv_type": svtype,
                "imprecise": bool(variant.INFO.get("IMPRECISE", False)),
                "length": float(length) if length and length > 0 else np.nan,
            }
            # cyvcf2 .genotypes: [allele_1, allele_2, ..., phased_bool]; -1 is a
            # missing allele. Re-serialize to the familiar GT string form.
            for sample_id, gt in zip(sample_ids, variant.genotypes):
                *alleles, phased = gt
                sep = "|" if phased else "/"
                row[sample_id] = sep.join("." if a < 0 else str(a) for a in alleles)
            rows.append(row)
        if n_skipped:
            log.warning("Skipped %d --vcf record(s) with unrecognized/missing SVTYPE", n_skipped)
        if not rows:
            log.warning("No usable SV records parsed from --vcf %s", vcf_path)
            return None
        log.info("Parsed %d SV record(s) across %d sample(s) from --vcf %s", len(rows), len(sample_ids), vcf_path)
        return pd.DataFrame(rows)
    except Exception as exc:  # noqa: BLE001 - best-effort parse, any failure -> caller falls back
        log.warning("Failed to parse --vcf %s: %s", vcf_path, exc)
        return None


def _find_column(columns: list[str], candidates: list[str]) -> str | None:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    return None


def load_real_haploblocks(path: Path) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:  # noqa: BLE001
        log.warning("Failed to parse haploblock table %s: %s", path, exc)
        return None
    chrom_col = _find_column(list(df.columns), ["chrom", "chr", "#chrom", "chromosome"])
    start_col = _find_column(list(df.columns), ["start", "chromstart", "block_start"])
    end_col = _find_column(list(df.columns), ["end", "chromend", "block_end"])
    if not (chrom_col and start_col and end_col):
        log.warning("Haploblock table %s is missing chrom/start/end columns", path)
        return None
    id_col = _find_column(list(df.columns), ["haploblock_id", "block_id", "id", "name"])
    cluster_col = _find_column(list(df.columns), ["n_clusters", "cluster_id", "cluster", "pop_cluster"])
    diff_col = _find_column(list(df.columns), ["cluster_diff_score", "diff_score", "fst"])
    if not cluster_col:
        log.warning("Haploblock table %s has no cluster-count column; population-cluster correlation (Stage 6) will have nothing to compare against", path)
    log.warning("Haploblock table %s has no SNP-density column; n_snps set to NaN -- Stage 5's block-architecture offset will need to handle missing values", path)
    out = pd.DataFrame(
        {
            "haploblock_id": df[id_col] if id_col else [f"hb{i}" for i in range(len(df))],
            "chrom": df[chrom_col],
            "start": df[start_col].astype(int),
            "end": df[end_col].astype(int),
            "n_snps": np.nan,
            "n_clusters": df[cluster_col] if cluster_col else np.nan,
            "hash_length": np.nan,
            "cluster_diff_score": df[diff_col] if diff_col else np.nan,
        }
    )
    return out


def load_real_sample_metadata(path: Path) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except Exception as exc:  # noqa: BLE001
        log.warning("Failed to parse sample metadata %s: %s", path, exc)
        return None
    sample_col = _find_column(list(df.columns), ["sample", "sample_id", "sampleid"])
    superpop_col = _find_column(list(df.columns), ["super_pop", "superpop", "superpopulation"])
    pop_col = _find_column(list(df.columns), ["pop", "population"])
    if not (sample_col and (superpop_col or pop_col)):
        log.warning("Sample metadata %s is missing sample/(super)population columns", path)
        return None
    superpop = df[superpop_col] if superpop_col else "UNKNOWN"
    # `population` is what Stages 4/6/7 group by; fall back to the superpopulation
    # label if the panel only carries the coarse grouping
    population = df[pop_col] if pop_col else superpop
    return pd.DataFrame(
        {"sample_id": df[sample_col], "population": population, "superpopulation": superpop}
    )


# --------------------------------------------------------------------------
# Chromosome-naming harmonization
# --------------------------------------------------------------------------

def normalize_chrom_column(df: pd.DataFrame) -> pd.DataFrame:
    """Force UCSC-style 'chrN' naming so tables from different sources can be joined.

    dbVar (and many NCBI/Ensembl-derived sources) use bare names ('1', 'X');
    data.haploblocks.org and this script's synthetic data use 'chr1'/'chrX'. Left
    unreconciled, an intersection between two such tables would silently match
    nothing on affected chromosomes despite both referring to the same physical
    sequence.
    """
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str).apply(lambda c: c if c.startswith("chr") else f"chr{c}")
    return df


# --------------------------------------------------------------------------
# Genome build harmonization
# --------------------------------------------------------------------------

def harmonize_build(
    df: pd.DataFrame, source_build: str, target_build: str, chain_file: str | None, label: str
) -> pd.DataFrame:
    if source_build == target_build:
        return df
    log.warning(
        "Genome build mismatch for %s: source is %s, pipeline target is %s",
        label,
        source_build,
        target_build,
    )
    if not chain_file:
        log.warning(
            "No --chain-file provided; proceeding WITHOUT liftover for %s -- "
            "its coordinates may not align with other inputs!",
            label,
        )
        return df
    try:
        from pyliftover import LiftOver

        lo = LiftOver(chain_file)
    except Exception as exc:  # noqa: BLE001
        log.warning("Could not load chain file %s (%s); skipping liftover for %s", chain_file, exc, label)
        return df

    lifted_starts, lifted_ends, n_failed = [], [], 0
    for _, row in df.iterrows():
        conv_start = lo.convert_coordinate(row["chrom"], int(row["start"]))
        conv_end = lo.convert_coordinate(row["chrom"], int(row["end"]))
        if conv_start and conv_end:
            lifted_starts.append(conv_start[0][1])
            lifted_ends.append(conv_end[0][1])
        else:
            lifted_starts.append(row["start"])
            lifted_ends.append(row["end"])
            n_failed += 1
    df = df.copy()
    df["start"], df["end"] = lifted_starts, lifted_ends
    if n_failed:
        log.warning("Liftover failed for %d/%d %s records; kept original coordinates for those", n_failed, len(df), label)
    else:
        log.info("Successfully lifted %d %s records to %s", len(df), label, target_build)
    return df


# --------------------------------------------------------------------------
# Synthetic fallback data (biologically plausible, not literally random)
# --------------------------------------------------------------------------

def generate_synthetic_haploblocks(n_blocks: int, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    base, remainder = divmod(n_blocks, len(CHROMS))
    block_idx = 0
    for chrom_i, chrom in enumerate(CHROMS):
        blocks_this_chrom = base + (1 if chrom_i < remainder else 0)
        # start each chromosome's covered span a few Mb in, so there is a
        # "before the first block" region for outside_block SVs to land in
        pos = int(rng.integers(1_000_000, 5_000_000))
        for _ in range(blocks_this_chrom):
            length = int(rng.integers(50_000, 500_000))
            start, end = pos, pos + length
            n_snps = int(max(10, rng.normal(length / 200, length / 1000)))
            n_clusters = int(rng.integers(2, 20))  # a block typically resolves into several haplotype clusters
            hash_length = int(rng.integers(10, 30))
            cluster_diff_score = float(rng.beta(2, 5))  # Fst-like: mostly low, occasionally high
            rows.append(
                {
                    "haploblock_id": f"hb{block_idx}",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "n_snps": n_snps,
                    "n_clusters": n_clusters,
                    "hash_length": hash_length,
                    "cluster_diff_score": round(cluster_diff_score, 4),
                }
            )
            pos = end  # blocks are contiguous within a chromosome, like real haploblocks.org data
            block_idx += 1
    return pd.DataFrame(rows)


def generate_synthetic_samples(n_samples: int, rng: np.random.Generator) -> pd.DataFrame:
    superpops = [SUPERPOPULATIONS[i % len(SUPERPOPULATIONS)] for i in range(n_samples)]
    # give each superpopulation two fake sub-populations so downstream code that
    # groups by `population` is genuinely exercised on >5 groups
    population = [f"{sp}{(i // len(SUPERPOPULATIONS)) % 2 + 1}" for i, sp in enumerate(superpops)]
    return pd.DataFrame(
        {
            "sample_id": [f"SAMP{i:04d}" for i in range(n_samples)],
            "population": population,
            "superpopulation": superpops,
        }
    )


def generate_synthetic_svs(
    n_svs: int,
    haploblocks: pd.DataFrame,
    samples: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    sample_pops = dict(zip(samples["sample_id"], samples["superpopulation"]))
    blocks_by_chrom = {
        c: g.sort_values("start").reset_index(drop=True)
        for c, g in haploblocks.groupby("chrom", sort=False)
    }
    chroms = list(blocks_by_chrom)
    rows = []
    for i in range(n_svs):
        sv_type = rng.choice(SV_TYPES, p=SV_TYPE_WEIGHTS)
        # interval_length places the SV on the reference; for INS this is just the
        # insertion point (1bp), matching how dbVar itself represents insertions --
        # the actual inserted-sequence length is a separate `sv_length` below
        interval_length = 1 if sv_type == "INS" else int(rng.lognormal(mean=8, sigma=1.2))
        interval_length = min(interval_length, 50_000)

        placement = rng.choice(["within", "boundary", "outside"], p=[0.6, 0.2, 0.2])
        chrom = chroms[int(rng.integers(0, len(chroms)))]
        blocks = blocks_by_chrom[chrom]

        if placement == "boundary" and len(blocks) >= 2:
            # straddle a shared edge: SV overlaps block k and block k+1
            k = int(rng.integers(0, len(blocks) - 1))
            edge = int(blocks.loc[k, "end"])  # == blocks.loc[k + 1, "start"] (contiguous)
            half = max(1, interval_length // 2)
            start = edge - half
            end = edge + max(1, interval_length - half)
        elif placement == "outside":
            if rng.random() < 0.5:  # before the first block on this chromosome
                end = int(blocks.loc[0, "start"]) - 1 - int(rng.integers(0, 20_000))
                start = end - interval_length
            else:  # after the last block
                start = int(blocks.iloc[-1]["end"]) + int(rng.integers(1, 20_000))
                end = start + interval_length
        else:  # within a single block, not touching either edge
            k = int(rng.integers(0, len(blocks)))
            b_start, b_end = int(blocks.loc[k, "start"]), int(blocks.loc[k, "end"])
            usable = (b_end - b_start) - interval_length
            start = b_start + 1 + int(rng.integers(0, max(1, usable - 1))) if usable > 2 else b_start
            end = min(start + interval_length, b_end)
        # plausible inserted-sequence length (Alu- to L1-scale), independent of
        # the 1bp reference interval -- mirrors real dbVar INS records, where
        # length lives with the call/insertion, not the reference coordinates
        sv_length = float(min(int(rng.lognormal(mean=6, sigma=1.0)), 20_000)) if sv_type == "INS" else float(interval_length)
        imprecise = bool(rng.random() < 0.28)  # matches the ~28% IMPRECISE fraction observed in real nstd152 data

        # allele-frequency profile: common / population-specific / rare
        sv_class = rng.choice(["common", "population_specific", "rare"], p=[0.5, 0.35, 0.15])
        af_by_pop = {}
        if sv_class == "common":
            base_af = rng.uniform(0.1, 0.4)
            af_by_pop = {p: max(0.01, base_af + rng.normal(0, 0.05)) for p in SUPERPOPULATIONS}
        elif sv_class == "population_specific":
            hot_pop = rng.choice(SUPERPOPULATIONS)
            af_by_pop = {p: rng.uniform(0.1, 0.3) if p == hot_pop else rng.uniform(0, 0.02) for p in SUPERPOPULATIONS}
        else:  # rare
            af_by_pop = {p: rng.uniform(0, 0.02) for p in SUPERPOPULATIONS}

        row = {
            "sv_id": f"sv{i}",
            "chrom": chrom,
            "start": start,
            "end": end,
            "sv_type": sv_type,
            "imprecise": imprecise,
            "length": sv_length,
        }
        for sample_id, pop in sample_pops.items():
            af = af_by_pop[pop]
            row[sample_id] = int(rng.binomial(2, af))
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def build_config(args, out_dir: Path, sources_used: dict) -> dict:
    return {
        "genome_build": args.target_genome_build,
        "data_sources": sources_used,
        "thresholds": {
            "af_common_threshold": args.af_common_threshold,
            "boundary_distance_bp": args.boundary_distance_bp,
            "min_sv_count_per_block": args.min_sv_count_per_block,
            "min_sv_length": args.min_sv_length,
            "max_sv_length": args.max_sv_length,
            "drop_imprecise": args.drop_imprecise,
        },
        "seeds": {
            "permutation_seed": args.seed,
            "umap_seed": args.seed,
        },
        "paths": {
            "sv_calls": str((out_dir / "sv_calls.tsv").resolve()),
            "haploblocks": str((out_dir / "haploblocks.tsv").resolve()),
            "sample_metadata": str((out_dir / "sample_metadata.tsv").resolve()),
        },
    }


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", default=None, help="Local path or URL to a standard multi-sample SV VCF (FORMAT/GT). Highest-priority SV source: builds sv_calls.tsv with the raw GT string per sample and derives sample_metadata.tsv from the VCF's sample list (population=UNKNOWN unless --panel-source is also given)")
    p.add_argument("--sv-source", default=None, help="Local path or URL to a single already-prepared SV VCF, genotypes collapsed to 0/1/2 dosage (overrides the dbVar auto-fetch below; --vcf overrides this)")
    p.add_argument("--sv-genome-build", default="GRCh38", help="Genome build of --sv-source, if given")
    p.add_argument("--dbvar-base-url", default=DBVAR_BASE_URL, help="Directory containing <study>.<build>.variant_{call,region}.vcf.gz(.tbi)")
    p.add_argument("--dbvar-study", default="nstd152")
    p.add_argument("--dbvar-build", default="GRCh38")
    p.add_argument("--skip-dbvar-download", action="store_true", help="Skip the automatic dbVar fetch; go straight to synthetic SV data (unless --sv-source is given)")
    p.add_argument("--haploblock-source", default=None, help="Local path or URL to a single already-prepared haploblock/cluster table (overrides the haploblocks.org auto-fetch below)")
    p.add_argument("--haploblock-genome-build", default="GRCh38", help="Genome build of --haploblock-source, if given")
    p.add_argument("--haploblock-hash-base-url", default=HAPLOBLOCK_HASH_BASE_URL)
    p.add_argument("--haploblock-chroms", default="all", help="Comma-separated chromosomes (e.g. chr21,chr22), or 'all' (default) for 1-22+X. One small aggregate-boundaries file is fetched per chromosome (23 requests total for 'all'), not one file per haploblock")
    p.add_argument("--skip-haploblock-download", action="store_true", help="Skip the automatic haploblocks.org fetch; go straight to synthetic haploblocks (unless --haploblock-source is given)")
    p.add_argument("--panel-source", default=None, help="Local path or URL to 1000 Genomes sample panel")
    p.add_argument("--target-genome-build", default="GRCh38")
    p.add_argument("--chain-file", default=None, help="UCSC-style chain file for liftover, if a build mismatch is found")
    p.add_argument("--out-dir", default="example_data", help="Output directory for tables + config.yaml")
    p.add_argument("--n-svs", type=int, default=300, help="Synthetic SV count (used if real data unavailable)")
    p.add_argument("--n-haploblocks", type=int, default=30, help="Synthetic haploblock count")
    p.add_argument("--n-samples", type=int, default=30, help="Synthetic sample count")
    p.add_argument("--boundary-distance-bp", type=int, default=5000, help="N: distance defining a boundary-crossing SV")
    p.add_argument("--af-common-threshold", type=float, default=0.05)
    p.add_argument("--min-sv-count-per-block", type=int, default=3)
    p.add_argument("--min-sv-length", type=int, default=50, help="Stage 1's size filter floor, bp (50bp is the conventional SV-vs-indel cutoff)")
    p.add_argument("--max-sv-length", type=int, default=5_000_000, help="Stage 1's size filter ceiling, bp")
    p.add_argument(
        "--drop-imprecise", action="store_true",
        help="Set thresholds.drop_imprecise: true in the emitted config so Stage 1 drops "
        "IMPRECISE-flagged calls. Off by default: callers flag nearly every inversion IMPRECISE, "
        "so dropping imprecise calls silently removes all INV events.",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    sources_used = {}

    # --- sample metadata (loaded first: SV synthesis needs sample -> population mapping) ---
    samples = None
    panel_path = fetch_or_load(args.panel_source, out_dir, "sample_metadata")
    if panel_path:
        samples = load_real_sample_metadata(panel_path)
        # no build harmonization here: a sample->superpopulation panel carries no genomic coordinates
    if samples is not None:
        log.info("Using real sample metadata (%d samples) from %s", len(samples), args.panel_source)
        sources_used["sample_metadata"] = "real"
    else:
        log.info("Using SYNTHETIC sample metadata (n=%d) -- no usable --panel-source given", args.n_samples)
        samples = generate_synthetic_samples(args.n_samples, rng)
        sources_used["sample_metadata"] = "synthetic"

    # --- haploblocks: manual --haploblock-source overrides the haploblocks.org auto-fetch ---
    haploblocks = None
    hb_path = fetch_or_load(args.haploblock_source, out_dir, "haploblocks")
    if hb_path:
        haploblocks = load_real_haploblocks(hb_path)
        if haploblocks is not None:
            haploblocks = harmonize_build(haploblocks, args.haploblock_genome_build, args.target_genome_build, args.chain_file, "haploblocks")
    elif not args.skip_haploblock_download:
        chroms = (
            ALL_1000G_CHROMS
            if args.haploblock_chroms.strip().lower() == "all"
            else [c.strip() for c in args.haploblock_chroms.split(",") if c.strip()]
        )
        haploblocks = fetch_haploblock_boundaries(args.haploblock_hash_base_url, chroms)
        if haploblocks is not None:
            # data.haploblocks.org is GRCh38-based; still route through the same harmonization
            # path as every other source for a consistent, non-silent build-mismatch check
            haploblocks = harmonize_build(haploblocks, "GRCh38", args.target_genome_build, args.chain_file, "haploblocks")
    if haploblocks is not None:
        log.info("Using real haploblock table (%d blocks) from %s", len(haploblocks), args.haploblock_source or args.haploblock_hash_base_url)
        sources_used["haploblocks"] = "real"
    else:
        log.info("Using SYNTHETIC haploblocks (n=%d) -- no usable haploblock source available", args.n_haploblocks)
        haploblocks = generate_synthetic_haploblocks(args.n_haploblocks, rng)
        sources_used["haploblocks"] = "synthetic"

    # --- SV calls: --vcf (raw GT) > --sv-source (dosage) > dbVar auto-fetch ---
    sv_calls = None
    sv_from_vcf = False
    if args.vcf:
        vcf_path = fetch_or_load(args.vcf, out_dir, "sv_calls_vcf")
        if vcf_path:
            sv_calls = load_vcf_sv_calls(vcf_path)
            if sv_calls is not None:
                sv_calls = harmonize_build(sv_calls, args.sv_genome_build, args.target_genome_build, args.chain_file, "sv_calls")
                sv_from_vcf = True
        if sv_calls is None:
            log.warning("--vcf %s yielded no usable calls; falling back to the other SV sources", args.vcf)

    sv_path = None if sv_from_vcf else fetch_or_load(args.sv_source, out_dir, "sv_calls")
    if sv_path:
        sv_calls = load_real_sv_calls(sv_path)
        if sv_calls is not None:
            sv_calls = harmonize_build(sv_calls, args.sv_genome_build, args.target_genome_build, args.chain_file, "sv_calls")
    elif not sv_from_vcf and not args.skip_dbvar_download:
        dbvar_paths = fetch_dbvar_study(args.dbvar_base_url, args.dbvar_study, args.dbvar_build, out_dir)
        if dbvar_paths:
            sv_calls = load_dbvar_sv_calls(dbvar_paths["call_vcf"], dbvar_paths["region_vcf"])
            if sv_calls is not None:
                sv_calls = harmonize_build(sv_calls, args.dbvar_build, args.target_genome_build, args.chain_file, "sv_calls")
    non_sample_cols = {"sv_id", "chrom", "start", "end", "sv_type", "imprecise", "length"}
    if sv_calls is not None:
        src_label = args.vcf if sv_from_vcf else (args.sv_source or f"{args.dbvar_base_url} ({args.dbvar_study})")
        log.info("Using real SV calls (%d records) from %s", len(sv_calls), src_label)
        sources_used["sv_calls"] = "real (--vcf)" if sv_from_vcf else "real"

        if sv_from_vcf and sources_used["sample_metadata"] != "real":
            # a plain VCF header gives sample ids but no population labels -- build
            # sample_metadata from the ids in the VCF so the sample_id columns agree
            # downstream; population is UNKNOWN until a --panel-source fills it in
            observed_samples = [c for c in sv_calls.columns if c not in non_sample_cols]
            samples = pd.DataFrame(
                {
                    "sample_id": observed_samples,
                    "population": "UNKNOWN",
                    "superpopulation": "UNKNOWN",
                }
            )
            sources_used["sample_metadata"] = "from --vcf sample ids (population=UNKNOWN)"
            log.warning(
                "--vcf gave %d sample(s) but no population labels; sample_metadata.tsv's "
                "population column is all UNKNOWN. Pass --panel-source to attach real "
                "population labels -- Stages 4/6/7 group by that column.",
                len(observed_samples),
            )
        elif sources_used["sample_metadata"] == "synthetic":
            # the synthetic sample panel's IDs (SAMPxxxx) won't match real dbVar sample
            # IDs (e.g. HG00512) -- rebuild sample_metadata from the IDs actually observed
            # in sv_calls so the two tables' sample_id columns agree downstream
            observed_samples = [c for c in sv_calls.columns if c not in non_sample_cols]
            unknown = [s for s in observed_samples if s not in HGSVC_TRIO_SUPERPOPULATIONS]
            if unknown:
                log.warning(
                    "%d observed SV sample id(s) have no known population mapping (not in the "
                    "hardcoded HGSVC trio table): %s -- labeled UNKNOWN", len(unknown), unknown,
                )
            samples = pd.DataFrame(
                {
                    "sample_id": observed_samples,
                    "population": [HGSVC_TRIO_POPULATIONS.get(s, "UNKNOWN") for s in observed_samples],
                    "superpopulation": [HGSVC_TRIO_SUPERPOPULATIONS.get(s, "UNKNOWN") for s in observed_samples],
                }
            )
            covered = sorted(set(samples["superpopulation"]) - {"UNKNOWN"})
            log.warning(
                "No --panel-source given: replaced the placeholder synthetic sample panel with %d "
                "sample(s) inferred from the real dbVar SV calls (hardcoded HGSVC trio lookup). "
                "This covers only %s -- population-specific classification cannot be evaluated for "
                "any superpopulation not listed. Pass --panel-source for the real IGSR panel to fix this.",
                len(samples), covered,
            )
            sources_used["sample_metadata"] = "real (inferred from dbVar SAMPLE ids)"
    else:
        log.info("Using SYNTHETIC SV calls (n=%d) -- no usable SV source available", args.n_svs)
        sv_calls = generate_synthetic_svs(args.n_svs, haploblocks, samples, rng)
        sources_used["sv_calls"] = "synthetic"

    sv_calls = normalize_chrom_column(sv_calls)
    haploblocks = normalize_chrom_column(haploblocks)

    sv_calls.to_csv(out_dir / "sv_calls.tsv", sep="\t", index=False)
    haploblocks.to_csv(out_dir / "haploblocks.tsv", sep="\t", index=False)
    samples.to_csv(out_dir / "sample_metadata.tsv", sep="\t", index=False)

    config = build_config(args, out_dir, sources_used)
    with open(out_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(config, fh, sort_keys=False)

    log.info("Wrote %d SVs, %d haploblocks, %d samples to %s", len(sv_calls), len(haploblocks), len(samples), out_dir.resolve())
    log.info("Data sources: %s", sources_used)
    log.info("Config written to %s", (out_dir / "config.yaml").resolve())


if __name__ == "__main__":
    main(sys.argv[1:])
