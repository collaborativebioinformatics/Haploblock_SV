"""Stage 5: per-haploblock SV-type enrichment.

Reads Stage 2's annotated sv_calls.tsv + haploblocks.tsv (paths from
--config) and asks, for every (haploblock, SV type) pair: does this block
carry significantly more (or fewer) SVs of this type than its length alone
would predict?

Method (deliberately simple -- negative-binomial / SNP-density / min-count
variants are parked in README.md "Descoped / future steps"):

  1. Assign each SV to its haploblock(s). `outside_block` SVs (no
     haploblock_id) are dropped; a `boundary_crossing` SV whose
     haploblock_id is comma-joined is counted once in EACH block it spans.
  2. observed[block, type]  = number of assigned SVs of that type.
  3. rate[type]             = (total assigned SVs of that type) /
                              (total length of all haploblocks)
  4. expected[block, type]  = rate[type] * length(block)
  5. p_value                = exact two-sided Poisson test of observed
                              against mean = expected
     (2 * min(P(X<=obs), P(X>=obs)), capped at 1).
  6. q_value                = Benjamini-Hochberg FDR across every
                              (block, type) cell tested.
  7. flagged                = q_value < --q-threshold (default 0.05).

The test grid is the full haploblocks x observed-types matrix, zeros
included -- so on a genome-wide run it is (n_haploblocks * n_types) rows and
most cells are observed=0. That is the honest denominator for "across all
haploblock x SV-type tests"; a min-count pre-filter (future work) would
trade some of that rigor for power.

Output: sv_type_enrichment.tsv with columns
  haploblock_id, sv_type, observed_count, expected_count, p_value, q_value, flagged
sorted by q_value ascending.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage5_type_enrichment")

try:
    from scipy.stats import false_discovery_control

    def bh_fdr(pvals: np.ndarray) -> np.ndarray:
        return false_discovery_control(np.asarray(pvals, dtype=float), method="bh")
except ImportError:  # older scipy
    from statsmodels.stats.multitest import multipletests

    def bh_fdr(pvals: np.ndarray) -> np.ndarray:
        return multipletests(np.asarray(pvals, dtype=float), method="fdr_bh")[1]


def load_config(config_path: Path) -> dict:
    with open(config_path) as fh:
        config = yaml.safe_load(fh)
    if "paths" not in config:
        raise ValueError(f"{config_path} is missing required top-level key 'paths'")
    for key in ("sv_calls", "haploblocks"):
        if key not in config["paths"]:
            raise ValueError(f"{config_path}: paths.{key} is required")
    return config


def _resolve(path_str: str, base_dir: Path) -> Path:
    p = Path(path_str)
    return p if p.is_absolute() else (base_dir / p)


def assign_svs_to_blocks(sv: pd.DataFrame) -> pd.DataFrame:
    """Return sv with one row per (SV, haploblock): drop rows with no
    haploblock_id, split comma-joined ids from boundary_crossing SVs."""
    hb_id = sv["haploblock_id"].fillna("").astype(str).str.strip()
    sv = sv.loc[hb_id != ""].copy()
    sv["haploblock_id"] = hb_id[hb_id != ""]
    sv = sv.assign(haploblock_id=sv["haploblock_id"].str.split(",")).explode("haploblock_id")
    sv["haploblock_id"] = sv["haploblock_id"].str.strip()
    return sv


def enrichment_table(sv: pd.DataFrame, hb: pd.DataFrame, q_threshold: float) -> pd.DataFrame:
    hb = hb.drop_duplicates("haploblock_id").copy()
    hb["length"] = hb["end"] - hb["start"]
    bad_len = hb["length"] <= 0
    if bad_len.any():
        log.warning("%d haploblock(s) have length <= 0 and are dropped: %s",
                    int(bad_len.sum()), hb.loc[bad_len, "haploblock_id"].tolist()[:10])
        hb = hb.loc[~bad_len]
    total_length = int(hb["length"].sum())

    assigned = assign_svs_to_blocks(sv)
    unknown = set(assigned["haploblock_id"]) - set(hb["haploblock_id"])
    if unknown:
        log.warning("%d assigned haploblock_id(s) not in the haploblock table, ignored: %s",
                    len(unknown), sorted(unknown)[:10])
        assigned = assigned[assigned["haploblock_id"].isin(hb["haploblock_id"])]

    sv_types = sorted(assigned["sv_type"].dropna().unique())
    if not sv_types:
        raise ValueError("no SVs could be assigned to a haploblock")

    observed = (
        assigned.groupby(["haploblock_id", "sv_type"]).size()
        .unstack("sv_type", fill_value=0)
        .reindex(index=hb["haploblock_id"], columns=sv_types, fill_value=0)
    )
    total_by_type = observed.sum(axis=0)
    rate_by_type = total_by_type / total_length
    log.info("overall SV rate per type (per bp, over %d bp of haploblock): %s",
             total_length, {t: round(float(r), 8) for t, r in rate_by_type.items()})

    length_by_block = hb.set_index("haploblock_id")["length"].reindex(observed.index)
    expected = np.outer(length_by_block.to_numpy(), rate_by_type.to_numpy())

    obs_flat = observed.to_numpy().ravel().astype(float)
    exp_flat = expected.ravel()
    # exact two-sided Poisson test of obs against mean = exp
    with np.errstate(invalid="ignore"):
        lower = stats.poisson.cdf(obs_flat, exp_flat)       # P(X <= obs)
        upper = stats.poisson.sf(obs_flat - 1, exp_flat)     # P(X >= obs)
    p_flat = np.minimum(1.0, 2.0 * np.minimum(lower, upper))
    p_flat = np.where(np.isfinite(p_flat), p_flat, 1.0)
    q_flat = bh_fdr(p_flat)

    rows = []
    block_ids = observed.index.to_numpy()
    n_types = len(sv_types)
    for i, block in enumerate(block_ids):
        for j, sv_type in enumerate(sv_types):
            k = i * n_types + j
            rows.append({
                "haploblock_id": block,
                "sv_type": sv_type,
                "observed_count": int(obs_flat[k]),
                "expected_count": float(exp_flat[k]),
                "p_value": float(p_flat[k]),
                "q_value": float(q_flat[k]),
                "flagged": bool(q_flat[k] < q_threshold),
            })
    return pd.DataFrame(rows).sort_values(["q_value", "haploblock_id", "sv_type"]).reset_index(drop=True)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="stage2_output/config.yaml", help="Stage 2's config.yaml (paths.sv_calls, paths.haploblocks)")
    p.add_argument("--out-dir", default="stage5_output", help="Output directory for the results table + this stage's own config.yaml")
    p.add_argument("--q-threshold", type=float, default=None, help="FDR q-value below which a (block, type) cell is flagged (default: config's thresholds.q_threshold, else 0.05)")
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)
    config_path = Path(args.config)
    config = load_config(config_path)
    base_dir = config_path.parent

    sv_path = _resolve(config["paths"]["sv_calls"], base_dir)
    hb_path = _resolve(config["paths"]["haploblocks"], base_dir)
    q_threshold = (
        args.q_threshold if args.q_threshold is not None
        else float(config.get("thresholds", {}).get("q_threshold", 0.05))
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sv = pd.read_csv(sv_path, sep="\t")
    hb = pd.read_csv(hb_path, sep="\t")
    for col, src in [("sv_type", sv_path), ("haploblock_id", sv_path)]:
        if col not in sv.columns:
            raise ValueError(f"{src} has no '{col}' column -- run Stage 2 first")

    result = enrichment_table(sv, hb, q_threshold)

    n_flag = int(result["flagged"].sum())
    log.info(
        "%d (haploblock x SV-type) cells tested (%d haploblocks x %d types); %d cell(s) with observed=0",
        len(result), result["haploblock_id"].nunique(), result["sv_type"].nunique(),
        int((result["observed_count"] == 0).sum()),
    )
    log.info("flagged (q < %.3g): %d", q_threshold, n_flag)
    if n_flag:
        show = result[result["flagged"]].head(15)
        for _, r in show.iterrows():
            log.info(
                "  %s  %s  observed=%d  expected=%.3f  q=%.2e",
                r["haploblock_id"], r["sv_type"], r["observed_count"], r["expected_count"], r["q_value"],
            )

    out_path = out_dir / "sv_type_enrichment.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    log.info("Wrote %d rows to %s", len(result), out_path.resolve())

    stage5_config = dict(config)
    stage5_config["paths"] = {
        "sv_calls": str(sv_path.resolve()),
        "haploblocks": str(hb_path.resolve()),
        "sv_type_enrichment": str(out_path.resolve()),
    }
    with open(out_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(stage5_config, fh, sort_keys=False)
    log.info("Config written to %s", (out_dir / "config.yaml").resolve())


if __name__ == "__main__":
    main(sys.argv[1:])
