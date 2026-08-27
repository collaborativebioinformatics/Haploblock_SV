"""Stage 4: common vs. population-specific SV classification.

Reads Stage 2's annotated sv_calls.tsv (paths from --config) and
sample_metadata.tsv, computes per-population alternate-allele frequency for
each SV, and labels the SV common / specific_to_population / other.

Populations come from sample_metadata.tsv's `population` column -- whatever
is in it. The five 1000G superpopulations are never hardcoded, and the
data.haploblocks.org cluster labels are deliberately NOT read here so that
Stage 6's SV-vs-cluster correlation stays non-circular.

Genotype columns in sv_calls.tsv may be either:
  - raw VCF GT strings ("0|1", "1/1", "./.")  -- Stage 0's --vcf path
  - 0/1/2 alt-allele dosage ints              -- Stage 0's dbVar/synthetic paths
Both are accepted. A missing genotype (`.` / `./.` / empty / NaN) is treated
as "not called", never as reference. Allele frequency for a population is
sum(alt alleles) / sum(called alleles) over that population's samples; a
population with no called genotype for an SV gets AF = NaN.

Output (tidy / long) -- one row per (SV x population):
  sv_id, sv_type, haploblock_id, chrom, start, end, position_class,
  population, n_samples, n_called,
  pop_has_data           True when >= --min-samples-per-pop of this
                          population's samples had a called genotype; only
                          populations with pop_has_data drive the category
  af                     alt-allele frequency in this population (NaN if
                          no called genotype)
  sv_category            common | specific_to_population | other
  specific_to_population  the population name when sv_category is
                          specific_to_population, else empty
  other_reason            absent_or_rare | one_pop_high_plus_intermediate_elsewhere
                          | insufficient_population_data   (empty otherwise)

Empty string cells (specific_to_population / other_reason for non-matching
rows) round-trip through pandas.read_csv as NaN -- a reader should
`.fillna("")` them. Any non-reference allele index in a GT string counts as
one alt allele (a fine approximation for biallelic SV callsets).

`sv_category` / `specific_to_population` / `other_reason` are constant across
the per-population rows of one SV.

Category rules (t = --af-threshold, default from config's
thresholds.af_common_threshold or 0.05; z = --absent-af-threshold, default
0.01; a population "has data" for an SV when at least --min-samples-per-pop
of its samples have a called genotype):
  - fewer than 2 populations have data           -> other (insufficient_population_data)
  - AF >= t in >= 2 populations                  -> common
  - AF >= t in exactly 1 population AND AF < z in
    every other population that has data          -> specific_to_population
  - otherwise (AF < t everywhere, or one high
    population with z <= AF < t elsewhere)        -> other

haploblock_id is passed through from Stage 2 unchanged: a single id for
within_block SVs, a comma-joined list for boundary_crossing SVs, empty for
outside_block SVs. Comma-joined ids are NOT split into separate rows.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage4_classify_af")

# columns of sv_calls.tsv that are NOT per-sample genotype columns
NON_GT_COLS = {
    "sv_id", "chrom", "start", "end", "sv_type", "imprecise", "length",
    "position_class", "haploblock_id",
}


def load_config(config_path: Path) -> dict:
    with open(config_path) as fh:
        config = yaml.safe_load(fh)
    if "paths" not in config:
        raise ValueError(f"{config_path} is missing required top-level key 'paths'")
    for key in ("sv_calls", "sample_metadata"):
        if key not in config["paths"]:
            raise ValueError(f"{config_path}: paths.{key} is required")
    return config


def _resolve(path_str: str, base_dir: Path) -> Path:
    """Resolve a config path relative to the config file's directory when it
    isn't absolute, so a checked-in example config can use bare filenames."""
    p = Path(path_str)
    return p if p.is_absolute() else (base_dir / p)


def genotype_columns(sv: pd.DataFrame) -> list[str]:
    return [c for c in sv.columns if c not in NON_GT_COLS]


def parse_genotypes(sv: pd.DataFrame, gt_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (alt, n_alleles): two [n_sv, n_sample] int frames -- the alt
    allele count and the called allele count for every SV x sample.

    A numeric column is read as a 0/1/2 alt dosage (called => 2 alleles).
    An object column is read as VCF GT strings; a cell that is just "0"/"1"/"2"
    (dosage encoded as text) is also accepted as a dosage.
    """
    alt = pd.DataFrame(0, index=sv.index, columns=gt_cols, dtype=int)
    n_alleles = pd.DataFrame(0, index=sv.index, columns=gt_cols, dtype=int)

    for c in gt_cols:
        col = sv[c]
        if pd.api.types.is_numeric_dtype(col):
            dosage = pd.to_numeric(col, errors="coerce")
            called = dosage.notna()
            alt[c] = dosage.where(called, 0).clip(0, 2).astype(int)
            n_alleles[c] = called.astype(int) * 2
            continue

        s = col.astype("string").str.strip()
        is_bare_dosage = s.str.fullmatch(r"[0-2]").fillna(False)
        bare_val = pd.to_numeric(s.where(is_bare_dosage), errors="coerce").fillna(0)
        alt_from_gt = s.str.count(r"[1-9]").fillna(0)          # any non-ref allele digit
        called_from_gt = s.str.count(r"[0-9]").fillna(0)        # any called allele digit
        alt[c] = np.where(is_bare_dosage, bare_val, alt_from_gt).astype(int)
        n_alleles[c] = np.where(is_bare_dosage, 2, called_from_gt).astype(int)

    return alt, n_alleles


def classify(
    sv: pd.DataFrame,
    meta: pd.DataFrame,
    gt_cols: list[str],
    af_threshold: float,
    absent_af_threshold: float,
    min_samples_per_pop: int,
) -> pd.DataFrame:
    """Return the tidy per-(SV x population) classification table."""
    alt, n_alleles = parse_genotypes(sv, gt_cols)

    pop_of = dict(zip(meta["sample_id"].astype(str), meta["population"].astype(str)))
    present = [s for s in gt_cols if s in pop_of]
    unknown = [s for s in gt_cols if s not in pop_of]
    if unknown:
        log.warning(
            "%d genotype column(s) absent from sample_metadata and ignored: %s",
            len(unknown), unknown[:10] + (["..."] if len(unknown) > 10 else []),
        )
    if not present:
        raise ValueError("no genotype column in sv_calls.tsv matches a sample_id in sample_metadata.tsv")

    pops = sorted({pop_of[s] for s in present})
    samples_by_pop = {p: [s for s in present if pop_of[s] == p] for p in pops}
    log.info(
        "%d population(s) from sample_metadata: %s",
        len(pops), {p: len(samples_by_pop[p]) for p in pops},
    )

    # per-population AF and called-sample count, one value per SV
    pop_af, pop_ncalled = {}, {}
    for p, samples in samples_by_pop.items():
        alt_sum = alt[samples].sum(axis=1).to_numpy()
        allele_sum = n_alleles[samples].sum(axis=1).to_numpy()
        with np.errstate(invalid="ignore", divide="ignore"):
            pop_af[p] = np.where(allele_sum > 0, alt_sum / allele_sum, np.nan)
        pop_ncalled[p] = (n_alleles[samples] > 0).sum(axis=1).to_numpy()

    categories, specifics, reasons = [], [], []
    for i in range(len(sv)):
        has_data = [p for p in pops if pop_ncalled[p][i] >= min_samples_per_pop]
        high = [p for p in has_data if pop_af[p][i] >= af_threshold]
        other_pops = [p for p in has_data if p not in high]

        if len(has_data) < 2:
            cat, spec, reason = "other", "", "insufficient_population_data"
        elif len(high) >= 2:
            cat, spec, reason = "common", "", ""
        elif len(high) == 1 and all(pop_af[p][i] < absent_af_threshold for p in other_pops):
            cat, spec, reason = "specific_to_population", high[0], ""
        elif len(high) == 1:
            cat, spec, reason = "other", "", "one_pop_high_plus_intermediate_elsewhere"
        else:
            cat, spec, reason = "other", "", "absent_or_rare"
        categories.append(cat)
        specifics.append(spec)
        reasons.append(reason)

    rows = []
    for i in range(len(sv)):
        r = sv.iloc[i]
        for p in pops:
            af = pop_af[p][i]
            rows.append({
                "sv_id": r["sv_id"],
                "sv_type": r["sv_type"],
                "haploblock_id": r.get("haploblock_id", ""),
                "chrom": r["chrom"],
                "start": r["start"],
                "end": r["end"],
                "position_class": r.get("position_class", ""),
                "population": p,
                "n_samples": len(samples_by_pop[p]),
                "n_called": int(pop_ncalled[p][i]),
                "pop_has_data": bool(pop_ncalled[p][i] >= min_samples_per_pop),
                "af": float(af) if np.isfinite(af) else np.nan,
                "sv_category": categories[i],
                "specific_to_population": specifics[i],
                "other_reason": reasons[i],
            })
    return pd.DataFrame(rows)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="stage2_output/config.yaml", help="Stage 2's config.yaml (paths.sv_calls, paths.sample_metadata)")
    p.add_argument("--out-dir", default="stage4_output", help="Output directory for the classification table + this stage's own config.yaml")
    p.add_argument("--af-threshold", type=float, default=None, help="AF at/above which a population 'has' the SV (default: config's thresholds.af_common_threshold, else 0.05)")
    p.add_argument("--absent-af-threshold", type=float, default=0.01, help="AF below which a population counts as 'near-zero' for the specific_to_population test (default 0.01)")
    p.add_argument("--min-samples-per-pop", type=int, default=2, help="Minimum called samples for a population to 'have data' for an SV (default 2)")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    config_path = Path(args.config)
    config = load_config(config_path)
    base_dir = config_path.parent

    sv_path = _resolve(config["paths"]["sv_calls"], base_dir)
    meta_path = _resolve(config["paths"]["sample_metadata"], base_dir)
    af_threshold = (
        args.af_threshold if args.af_threshold is not None
        else float(config.get("thresholds", {}).get("af_common_threshold", 0.05))
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sv = pd.read_csv(sv_path, sep="\t")
    meta = pd.read_csv(meta_path, sep="\t")
    if "population" not in meta.columns:
        raise ValueError(f"{meta_path} has no 'population' column -- Stages 4/6/7 group by it")
    if "haploblock_id" in sv.columns:
        sv["haploblock_id"] = sv["haploblock_id"].fillna("")

    gt_cols = genotype_columns(sv)
    if not gt_cols:
        raise ValueError(f"{sv_path} has no per-sample genotype columns")

    result = classify(
        sv, meta, gt_cols,
        af_threshold=af_threshold,
        absent_af_threshold=args.absent_af_threshold,
        min_samples_per_pop=args.min_samples_per_pop,
    )

    per_sv = result.drop_duplicates("sv_id")
    cat_counts = per_sv["sv_category"].value_counts()
    log.info(
        "Classified %d SV(s) (t=%.3g, near-zero<%.3g, min %d called samples/pop): %s",
        len(per_sv), af_threshold, args.absent_af_threshold, args.min_samples_per_pop,
        dict(cat_counts),
    )
    n_other = int(cat_counts.get("other", 0))
    if len(per_sv):
        log.info("'other' fraction: %.2f", n_other / len(per_sv))
    if n_other:
        log.info("'other' breakdown: %s", dict(per_sv.loc[per_sv["sv_category"] == "other", "other_reason"].value_counts()))

    out_path = out_dir / "sv_af_classification.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    log.info(
        "Wrote %d rows (%d SVs x %d populations) to %s",
        len(result), len(per_sv), result["population"].nunique(), out_path.resolve(),
    )

    stage4_config = dict(config)
    stage4_config["paths"] = {
        "sv_calls": str(sv_path.resolve()),
        "sample_metadata": str(meta_path.resolve()),
        "sv_af_classification": str(out_path.resolve()),
    }
    if "haploblocks" in config["paths"]:
        stage4_config["paths"]["haploblocks"] = str(_resolve(config["paths"]["haploblocks"], base_dir).resolve())
    with open(out_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(stage4_config, fh, sort_keys=False)
    log.info("Config written to %s", (out_dir / "config.yaml").resolve())


if __name__ == "__main__":
    main(sys.argv[1:])
