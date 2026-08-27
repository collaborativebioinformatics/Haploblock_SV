"""Tests for pipeline/stage0_ingest.py's synthetic-data (offline) path.

Runs the script as a subprocess with --skip-dbvar-download and
--skip-haploblock-download so nothing here touches the network, then checks
its outputs against the config.yaml/table contract every later stage relies on.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

STAGE0_SCRIPT = Path(__file__).resolve().parent.parent / "pipeline" / "stage0_ingest.py"
SUPERPOPULATIONS = {"AFR", "AMR", "EAS", "EUR", "SAS"}
SV_TYPES = {"DEL", "DUP", "INV", "INS"}


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("stage0_synthetic")
    result = subprocess.run(
        [
            sys.executable, str(STAGE0_SCRIPT),
            "--skip-dbvar-download",
            "--skip-haploblock-download",
            "--out-dir", str(out_dir),
            "--n-svs", "50",
            "--n-haploblocks", "6",
            "--n-samples", "10",
            "--seed", "1",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, f"stage0_ingest.py failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    return out_dir


def test_config_yaml_has_expected_keys(synthetic_run):
    config_path = synthetic_run / "config.yaml"
    assert config_path.exists()
    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    assert config["genome_build"] == "GRCh38"
    assert set(config["thresholds"]) == {
        "af_common_threshold", "boundary_distance_bp", "min_sv_count_per_block",
        "min_sv_length", "max_sv_length", "drop_imprecise",
    }
    assert set(config["seeds"]) == {"permutation_seed", "umap_seed"}
    assert set(config["paths"]) == {"sv_calls", "haploblocks", "sample_metadata"}
    assert config["data_sources"] == {
        "sample_metadata": "synthetic",
        "haploblocks": "synthetic",
        "sv_calls": "synthetic",
    }
    for path_str in config["paths"].values():
        assert Path(path_str).exists(), f"config.yaml references a missing file: {path_str}"


def test_sv_calls_table_shape_and_columns(synthetic_run):
    sv = pd.read_csv(synthetic_run / "sv_calls.tsv", sep="\t")
    assert isinstance(sv, pd.DataFrame)
    assert len(sv) > 0
    for col in ["sv_id", "chrom", "start", "end", "sv_type", "imprecise", "length"]:
        assert col in sv.columns
    assert set(sv["sv_type"].unique()).issubset(SV_TYPES)
    assert (sv["end"] >= sv["start"]).all()
    assert sv["imprecise"].dtype == bool
    assert (sv["length"].dropna() > 0).all()


def test_haploblocks_table_shape_and_columns(synthetic_run):
    hb = pd.read_csv(synthetic_run / "haploblocks.tsv", sep="\t")
    assert isinstance(hb, pd.DataFrame)
    assert len(hb) > 0
    for col in ["haploblock_id", "chrom", "start", "end", "n_snps", "n_clusters", "cluster_diff_score"]:
        assert col in hb.columns
    assert (hb["end"] > hb["start"]).all()


def test_sample_metadata_superpopulations(synthetic_run):
    samples = pd.read_csv(synthetic_run / "sample_metadata.tsv", sep="\t")
    assert len(samples) > 0
    assert "superpopulation" in samples.columns
    assert set(samples["superpopulation"].unique()).issubset(SUPERPOPULATIONS)
    # `population` is the fine-grained label Stages 4/6/7 group by
    assert "population" in samples.columns
    assert samples["population"].notna().all()


# --------------------------------------------------------------------------
# --vcf path: build the tables from a standard multi-sample VCF
# --------------------------------------------------------------------------

EXAMPLE_VCF = STAGE0_SCRIPT.parent.parent / "example_data" / "example_cohort.vcf"


@pytest.fixture(scope="module")
def vcf_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("stage0_vcf")
    result = subprocess.run(
        [
            sys.executable, str(STAGE0_SCRIPT),
            "--vcf", str(EXAMPLE_VCF),
            "--skip-dbvar-download", "--skip-haploblock-download",
            "--out-dir", str(out_dir),
        ],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, f"stage0 --vcf failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    return out_dir


def test_vcf_sv_calls_keep_raw_genotype_strings(vcf_run):
    sv = pd.read_csv(vcf_run / "sv_calls.tsv", sep="\t")
    sample_cols = [c for c in sv.columns if c not in
                   {"sv_id", "chrom", "start", "end", "sv_type", "imprecise", "length"}]
    assert sample_cols == ["SAMP_A", "SAMP_B", "SAMP_C"]
    # every per-sample cell is a GT string ("0/1", "1|1", "./."), not a bare dosage int
    gt_cells = sv[sample_cols].astype(str).to_numpy().ravel()
    assert all(("/" in c) or ("|" in c) for c in gt_cells)
    assert {"|", "."}.issubset(set("".join(gt_cells)))  # phasing + missingness both survived


def test_vcf_keeps_inversion_and_defaults_to_keeping_imprecise(vcf_run):
    sv = pd.read_csv(vcf_run / "sv_calls.tsv", sep="\t")
    inv = sv[sv["sv_type"] == "INV"]
    assert len(inv) == 1 and bool(inv["imprecise"].iloc[0]) is True
    with open(vcf_run / "config.yaml") as fh:
        config = yaml.safe_load(fh)
    assert config["thresholds"]["drop_imprecise"] is False


def test_vcf_sample_metadata_from_header(vcf_run):
    samples = pd.read_csv(vcf_run / "sample_metadata.tsv", sep="\t")
    assert list(samples["sample_id"]) == ["SAMP_A", "SAMP_B", "SAMP_C"]
    assert set(samples["population"]) == {"UNKNOWN"}
