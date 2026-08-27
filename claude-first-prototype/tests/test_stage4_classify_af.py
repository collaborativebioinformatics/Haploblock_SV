"""Tests for pipeline/stage4_classify_af.py.

  - end-to-end on the checked-in example_data/stage4_example/ (GT strings,
    4 populations, one SV per target category)
  - a direct-call check that 0/1/2 dosage-int genotype columns are handled
    the same way as GT strings.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
STAGE4_SCRIPT = REPO_ROOT / "pipeline" / "stage4_classify_af.py"
EXAMPLE_DIR = REPO_ROOT / "example_data" / "stage4_example"

# import the module for direct-call tests
_spec = importlib.util.spec_from_file_location("stage4_classify_af", STAGE4_SCRIPT)
stage4 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stage4)


@pytest.fixture(scope="module")
def example_output(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("stage4_out")
    result = subprocess.run(
        [sys.executable, str(STAGE4_SCRIPT),
         "--config", str(EXAMPLE_DIR / "config.yaml"), "--out-dir", str(out_dir)],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, f"stage4 failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    df = pd.read_csv(out_dir / "sv_af_classification.tsv", sep="\t")
    return df


def _cat(df, sv_id):
    return df.loc[df["sv_id"] == sv_id, "sv_category"].iloc[0]


def test_row_count_is_svs_times_populations(example_output):
    assert len(example_output) == 5 * 4
    assert set(example_output["population"]) == {"popA", "popB", "popC", "popD"}


def test_common_sv(example_output):
    assert _cat(example_output, "sv_common") == "common"
    af = example_output.set_index(["sv_id", "population"])["af"]
    assert af[("sv_common", "popA")] == pytest.approx(0.3)
    assert af[("sv_common", "popB")] == pytest.approx(0.3)
    assert af[("sv_common", "popC")] == pytest.approx(0.1)


def test_population_specific_svs_record_the_population(example_output):
    b = example_output[example_output["sv_id"] == "sv_specific_popB"].iloc[0]
    assert b["sv_category"] == "specific_to_population"
    assert b["specific_to_population"] == "popB"
    a = example_output[example_output["sv_id"] == "sv_specific_popA"].iloc[0]
    assert a["sv_category"] == "specific_to_population"
    assert a["specific_to_population"] == "popA"
    # popB private SV: AF is high only in popB
    afB = example_output.set_index(["sv_id", "population"])["af"][("sv_specific_popB", "popB")]
    assert afB == pytest.approx(0.5)


def test_other_categories_carry_a_reason(example_output):
    rare = example_output[example_output["sv_id"] == "sv_rare"].iloc[0]
    assert rare["sv_category"] == "other"
    assert rare["other_reason"] == "absent_or_rare"
    low = example_output[example_output["sv_id"] == "sv_lowdata"].iloc[0]
    assert low["sv_category"] == "other"
    assert low["other_reason"] == "insufficient_population_data"


def test_underpowered_population_is_flagged_not_dropped(example_output):
    # popD has a single sample: it appears in the table but never "has data"
    popd = example_output[example_output["population"] == "popD"]
    assert len(popd) == 5
    assert not popd["pop_has_data"].any()


def test_missing_genotypes_are_not_counted_as_reference(example_output):
    # sv_rare has one './.' in popA -> 4 of 5 samples called, AF still 0
    row = example_output[(example_output["sv_id"] == "sv_rare") & (example_output["population"] == "popA")].iloc[0]
    assert row["n_called"] == 4
    assert row["af"] == pytest.approx(0.0)


def test_dosage_int_genotypes_match_gt_strings():
    """A 0/1/2 dosage matrix must classify the same as the equivalent GT strings."""
    meta = pd.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4", "s5", "s6"],
        "population": ["P1", "P1", "P1", "P2", "P2", "P2"],
    })
    base = {"sv_id": "x", "sv_type": "DEL", "chrom": "chr1", "start": 1, "end": 2,
            "haploblock_id": "hb", "position_class": "within_block"}
    dosage = pd.DataFrame([{**base, "s1": 1, "s2": 1, "s3": 0, "s4": 0, "s5": 0, "s6": 0}])
    gt = pd.DataFrame([{**base, "s1": "0/1", "s2": "0|1", "s3": "0/0", "s4": "0/0", "s5": "0/0", "s6": "./."}])
    gt_cols = ["s1", "s2", "s3", "s4", "s5", "s6"]

    out_d = stage4.classify(dosage, meta, gt_cols, af_threshold=0.05, absent_af_threshold=0.01, min_samples_per_pop=2)
    out_g = stage4.classify(gt, meta, gt_cols, af_threshold=0.05, absent_af_threshold=0.01, min_samples_per_pop=2)

    afd = out_d.set_index("population")["af"]
    assert afd["P1"] == pytest.approx(2 / 6)   # 2 alt alleles / 6 called
    assert afd["P2"] == pytest.approx(0.0)
    # GT frame: P2 has one './.' so only 4 called alleles there, still AF 0
    afg = out_g.set_index("population")["af"]
    assert afg["P1"] == pytest.approx(2 / 6)
    assert afg["P2"] == pytest.approx(0.0)
    assert out_d["sv_category"].iloc[0] == out_g["sv_category"].iloc[0] == "specific_to_population"
    assert out_d["specific_to_population"].iloc[0] == "P1"
