"""Tests for pipeline/stage5_type_enrichment.py.

  - direct-call unit checks (SV->block assignment; a deterministic spike vs.
    length-proportional counts)
  - end-to-end on the checked-in example_data/stage5_example/
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
STAGE5_SCRIPT = REPO_ROOT / "pipeline" / "stage5_type_enrichment.py"
EXAMPLE_DIR = REPO_ROOT / "example_data" / "stage5_example"

_spec = importlib.util.spec_from_file_location("stage5_type_enrichment", STAGE5_SCRIPT)
stage5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stage5)


# --------------------------------------------------------------------------
# unit: SV -> haploblock assignment
# --------------------------------------------------------------------------

def test_assign_drops_outside_and_explodes_boundary_crossing():
    sv = pd.DataFrame([
        {"sv_id": "a", "sv_type": "DEL", "haploblock_id": "hb1", "position_class": "within_block"},
        {"sv_id": "b", "sv_type": "INV", "haploblock_id": "hb1,hb2", "position_class": "boundary_crossing"},
        {"sv_id": "c", "sv_type": "INS", "haploblock_id": "", "position_class": "outside_block"},
        {"sv_id": "d", "sv_type": "INS", "haploblock_id": np.nan, "position_class": "outside_block"},
    ])
    out = stage5.assign_svs_to_blocks(sv)
    assert list(out["sv_id"]) == ["a", "b", "b"]                 # c, d dropped; b in two blocks
    assert set(out.loc[out["sv_id"] == "b", "haploblock_id"]) == {"hb1", "hb2"}


# --------------------------------------------------------------------------
# unit: a deterministic spike is flagged, proportional counts are not
# --------------------------------------------------------------------------

def test_enrichment_flags_only_the_spike():
    # b_big is large enough to anchor the background rate; counts elsewhere are
    # length-proportional to it, except a small INV spike in the tiny block.
    hb = pd.DataFrame([
        {"haploblock_id": "b_small", "chrom": "chr1", "start": 0, "end": 1_000},
        {"haploblock_id": "b_mid", "chrom": "chr1", "start": 1_000, "end": 21_000},
        {"haploblock_id": "b_big", "chrom": "chr1", "start": 21_000, "end": 221_000},
    ])
    rows = []

    def block_svs(block_id, n_del, n_inv):
        for _ in range(n_del):
            rows.append({"sv_id": f"s{len(rows)}", "sv_type": "DEL", "haploblock_id": block_id, "position_class": "within_block"})
        for _ in range(n_inv):
            rows.append({"sv_id": f"s{len(rows)}", "sv_type": "INV", "haploblock_id": block_id, "position_class": "within_block"})

    # DEL ~ 1 / 500 bp, INV ~ 1 / 1000 bp, proportional to each block's length
    block_svs("b_small", n_del=2, n_inv=0)
    block_svs("b_mid", n_del=40, n_inv=20)
    block_svs("b_big", n_del=400, n_inv=200)
    # artificial INV spike in the tiny block (expected ~1, observed 8)
    block_svs("b_small", n_del=0, n_inv=8)

    result = stage5.enrichment_table(pd.DataFrame(rows), hb, q_threshold=0.05)

    flagged = result[result["flagged"]]
    assert list(flagged[["haploblock_id", "sv_type"]].itertuples(index=False, name=None)) == [("b_small", "INV")]
    # per type, sum(expected) == sum(observed) (rate is total / total length)
    agg = result.groupby("sv_type")[["observed_count", "expected_count"]].sum()
    assert np.allclose(agg["observed_count"], agg["expected_count"])


def test_length_proportional_draws_not_flagged_only_injected_spike_is():
    """Several haploblocks of varying length with SV counts drawn Poisson-
    proportional to length (fixed seed) -> none flagged after the length
    adjustment; one block gets an artificial INV inflation -> flagged after FDR."""
    rng = np.random.default_rng(0)
    blocks = [
        ("hb_a", 5_000), ("hb_b", 10_000), ("hb_c", 15_000), ("hb_d", 25_000),
        ("hb_e", 40_000), ("hb_f", 60_000), ("hb_g", 90_000), ("hb_anchor", 250_000),
    ]
    hb_rows, pos = [], 0
    for bid, length in blocks:
        hb_rows.append({"haploblock_id": bid, "chrom": "chr1", "start": pos, "end": pos + length})
        pos += length
    hb = pd.DataFrame(hb_rows)

    rates = {"DEL": 1 / 500, "INS": 1 / 1_200, "INV": 1 / 1_500}   # per bp
    sv_rows = []
    for bid, length in blocks:
        for sv_type, rate in rates.items():
            for _ in range(rng.poisson(rate * length)):
                sv_rows.append({"sv_id": f"s{len(sv_rows)}", "sv_type": sv_type,
                                "haploblock_id": bid, "position_class": "within_block"})
    # artificial INV inflation in one mid-size block
    for _ in range(15):
        sv_rows.append({"sv_id": f"s{len(sv_rows)}", "sv_type": "INV",
                        "haploblock_id": "hb_c", "position_class": "within_block"})

    result = stage5.enrichment_table(pd.DataFrame(sv_rows), hb, q_threshold=0.05)

    flagged = set(result.loc[result["flagged"], ["haploblock_id", "sv_type"]]
                  .itertuples(index=False, name=None))
    assert flagged == {("hb_c", "INV")}, flagged

    spike = result[(result["haploblock_id"] == "hb_c") & (result["sv_type"] == "INV")].iloc[0]
    assert spike["observed_count"] > spike["expected_count"]     # over-, not under-enriched
    assert spike["q_value"] < 0.05
    # every length-proportional (block x type) cell stays well above the cutoff
    proportional = result[~((result["haploblock_id"] == "hb_c") & (result["sv_type"] == "INV"))]
    assert (proportional["q_value"] > 0.05).all()


# --------------------------------------------------------------------------
# end-to-end on the bundled example
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def example_output(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("stage5_out")
    r = subprocess.run(
        [sys.executable, str(STAGE5_SCRIPT),
         "--config", str(EXAMPLE_DIR / "config.yaml"), "--out-dir", str(out_dir)],
        capture_output=True, text=True, timeout=60,
    )
    assert r.returncode == 0, f"stage5 failed:\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    return pd.read_csv(out_dir / "sv_type_enrichment.tsv", sep="\t")


def test_example_shape_and_columns(example_output):
    assert list(example_output.columns) == [
        "haploblock_id", "sv_type", "observed_count", "expected_count",
        "p_value", "q_value", "flagged",
    ]
    assert len(example_output) == 6 * 3  # 6 haploblocks x {DEL, INS, INV}


def test_example_flags_only_the_inv_spike(example_output):
    flagged = example_output[example_output["flagged"]]
    assert len(flagged) == 1
    row = flagged.iloc[0]
    assert (row["haploblock_id"], row["sv_type"]) == ("hb6_small", "INV")
    assert row["observed_count"] == 6
    assert row["expected_count"] < 0.5
    assert row["q_value"] < 0.05


def test_example_length_proportional_cells_not_flagged(example_output):
    prop = example_output[~(
        (example_output["haploblock_id"] == "hb6_small") & (example_output["sv_type"] == "INV")
    )]
    assert not prop["flagged"].any()
