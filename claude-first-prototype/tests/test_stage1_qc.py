"""Tests for pipeline/stage1_qc.py.

Chains off Stage 0's synthetic (offline) output rather than touching the
network, then checks the filter/dedup/validation behavior described in
stage1_qc.py's module docstring.
"""

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
STAGE0_SCRIPT = REPO_ROOT / "pipeline" / "stage0_ingest.py"
STAGE1_SCRIPT = REPO_ROOT / "pipeline" / "stage1_qc.py"


def run(script, args, timeout=60):
    result = subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True, timeout=timeout)
    assert result.returncode == 0, f"{script.name} failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    return result


@pytest.fixture(scope="module")
def stage0_synthetic(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("stage0_synthetic")
    run(
        STAGE0_SCRIPT,
        [
            "--skip-dbvar-download", "--skip-haploblock-download",
            "--out-dir", str(out_dir),
            "--n-svs", "200", "--n-haploblocks", "8", "--n-samples", "10",
            "--seed", "7",
        ],
    )
    return out_dir


@pytest.fixture(scope="module")
def stage1_output(stage0_synthetic, tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("stage1_output")
    run(STAGE1_SCRIPT, ["--config", str(stage0_synthetic / "config.yaml"), "--out-dir", str(out_dir)])
    return out_dir


def test_qc_report_counts_are_consistent(stage1_output, stage0_synthetic):
    with open(stage1_output / "qc_report.json") as fh:
        report = json.load(fh)

    raw = pd.read_csv(stage0_synthetic / "sv_calls.tsv", sep="\t")
    assert report["sv_calls"]["input_rows"] == len(raw)

    s = report["sv_calls"]
    assert s["input_rows"] - s["dropped_imprecise"] == s["rows_after_confidence_filter"]
    assert s["rows_after_confidence_filter"] - s["dropped_out_of_size_range"] == s["rows_after_size_filter"]
    assert s["rows_after_size_filter"] - s["exact_duplicates_removed"] == s["output_rows"]

    cleaned = pd.read_csv(stage1_output / "sv_calls.tsv", sep="\t")
    assert len(cleaned) == s["output_rows"]
    assert report["haploblocks"]["validation"].startswith("passed")


def test_cleaned_sv_calls_satisfy_filters(stage1_output, stage0_synthetic):
    with open(stage0_synthetic / "config.yaml") as fh:
        config = yaml.safe_load(fh)
    thresholds = config["thresholds"]

    sv = pd.read_csv(stage1_output / "sv_calls.tsv", sep="\t")
    assert len(sv) > 0
    assert not sv["imprecise"].any()
    in_range = sv["length"].isna() | sv["length"].between(thresholds["min_sv_length"], thresholds["max_sv_length"])
    assert in_range.all()
    assert (sv["start"] <= sv["end"]).all()
    assert not sv.duplicated(subset=["chrom", "start", "end", "sv_type", "length"]).any()


def test_haploblocks_passed_through_unchanged(stage1_output, stage0_synthetic):
    raw_hb = pd.read_csv(stage0_synthetic / "haploblocks.tsv", sep="\t")
    cleaned_hb = pd.read_csv(stage1_output / "haploblocks.tsv", sep="\t")
    pd.testing.assert_frame_equal(raw_hb, cleaned_hb)


def test_stage1_config_points_at_its_own_output(stage1_output):
    with open(stage1_output / "config.yaml") as fh:
        config = yaml.safe_load(fh)
    for key in ("sv_calls", "haploblocks", "sample_metadata"):
        path = Path(config["paths"][key])
        assert path.parent == stage1_output
        assert path.exists()


def test_overlapping_haploblocks_raise_with_details(stage0_synthetic, tmp_path):
    bad_dir = tmp_path / "bad_input"
    bad_dir.mkdir()
    hb = pd.read_csv(stage0_synthetic / "haploblocks.tsv", sep="\t")
    same_chrom_idx = hb[hb["chrom"] == hb["chrom"].iloc[0]].index[:2]
    assert len(same_chrom_idx) == 2, "fixture needs >=2 blocks on one chromosome for this test"
    hb.loc[same_chrom_idx[1], "start"] = hb.loc[same_chrom_idx[0], "end"] - 100
    hb.to_csv(bad_dir / "haploblocks.tsv", sep="\t", index=False)

    with open(stage0_synthetic / "config.yaml") as fh:
        config = yaml.safe_load(fh)
    config["paths"]["haploblocks"] = str(bad_dir / "haploblocks.tsv")
    with open(bad_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(config, fh)

    result = subprocess.run(
        [sys.executable, str(STAGE1_SCRIPT), "--config", str(bad_dir / "config.yaml"), "--out-dir", str(tmp_path / "out")],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode != 0
    assert "Haploblock BED validation failed" in result.stderr
    assert "overlaps" in result.stderr
