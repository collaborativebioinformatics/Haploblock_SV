"""Tests for pipeline/stage2_intersect.py.

Two kinds of coverage:
  - a hand-crafted tiny table (two adjacent blocks + a lone block on another
    chromosome) that exercises within_block / boundary_crossing (both the
    "near one edge" and "spans two blocks" cases) / outside_block exactly,
    per the classification rules in stage2_intersect.py's module docstring.
  - an integration check chaining off Stage 0 -> Stage 1's synthetic
    (offline) output, to confirm the --config wiring actually works.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
STAGE0_SCRIPT = REPO_ROOT / "pipeline" / "stage0_ingest.py"
STAGE1_SCRIPT = REPO_ROOT / "pipeline" / "stage1_qc.py"
STAGE2_SCRIPT = REPO_ROOT / "pipeline" / "stage2_intersect.py"


def run(script, args, timeout=60):
    result = subprocess.run([sys.executable, str(script), *args], capture_output=True, text=True, timeout=timeout)
    assert result.returncode == 0, f"{script.name} failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    return result


# --------------------------------------------------------------------------
# Hand-crafted classification test
# --------------------------------------------------------------------------

BOUNDARY_BP = 100


@pytest.fixture(scope="module")
def handcrafted_output(tmp_path_factory):
    in_dir = tmp_path_factory.mktemp("stage2_handcrafted_in")

    haploblocks = pd.DataFrame(
        [
            {"haploblock_id": "hbA", "chrom": "chr1", "start": 1000, "end": 2000, "n_snps": 100, "n_clusters": 3, "hash_length": 20, "cluster_diff_score": 0.1},
            {"haploblock_id": "hbB", "chrom": "chr1", "start": 2000, "end": 3000, "n_snps": 100, "n_clusters": 3, "hash_length": 20, "cluster_diff_score": 0.1},
            # isolated: >= BOUNDARY_BP away from both hbA and hbB, to test a genuinely
            # one-sided "near edge" case distinct from the hbA/hbB shared (contiguous) edge
            {"haploblock_id": "hbD", "chrom": "chr1", "start": 3500, "end": 4000, "n_snps": 100, "n_clusters": 3, "hash_length": 20, "cluster_diff_score": 0.1},
            {"haploblock_id": "hbC", "chrom": "chr2", "start": 5000, "end": 6000, "n_snps": 100, "n_clusters": 3, "hash_length": 20, "cluster_diff_score": 0.1},
        ]
    )
    sv_calls = pd.DataFrame(
        [
            # fully inside hbA, far (>=100bp) from both edges -> within_block
            {"sv_id": "sv_within", "chrom": "chr1", "start": 1400, "end": 1600, "sv_type": "DEL", "imprecise": False, "length": 200, "SAMP0": 1},
            # straddles the hbA/hbB shared edge at 2000 -> boundary_crossing, both ids
            {"sv_id": "sv_span", "chrom": "chr1", "start": 1990, "end": 2010, "sv_type": "DEL", "imprecise": False, "length": 20, "SAMP0": 1},
            # fully inside hbA, 40bp from its right edge (< BOUNDARY_BP) -- hbB starts
            # exactly there (contiguous), so hbB is also "near" -> boundary_crossing, both ids
            {"sv_id": "sv_near_shared_edge", "chrom": "chr1", "start": 1900, "end": 1960, "sv_type": "DUP", "imprecise": False, "length": 60, "SAMP0": 0},
            # fully inside hbD, 40bp from its left edge, with no neighboring block within
            # BOUNDARY_BP on either side -> boundary_crossing, hbD only
            {"sv_id": "sv_near_isolated_edge", "chrom": "chr1", "start": 3540, "end": 3560, "sv_type": "DUP", "imprecise": False, "length": 20, "SAMP0": 0},
            # far from every block on chr1 -> outside_block
            {"sv_id": "sv_far", "chrom": "chr1", "start": 100, "end": 200, "sv_type": "INV", "imprecise": False, "length": 100, "SAMP0": 0},
            # chromosome with zero haploblocks -> outside_block
            {"sv_id": "sv_no_chrom", "chrom": "chr3", "start": 1, "end": 10, "sv_type": "INS", "imprecise": False, "length": 9, "SAMP0": 1},
        ]
    )
    sample_metadata = pd.DataFrame([{"sample_id": "SAMP0", "superpopulation": "EUR"}])

    haploblocks.to_csv(in_dir / "haploblocks.tsv", sep="\t", index=False)
    sv_calls.to_csv(in_dir / "sv_calls.tsv", sep="\t", index=False)
    sample_metadata.to_csv(in_dir / "sample_metadata.tsv", sep="\t", index=False)
    config = {
        "genome_build": "GRCh38",
        "data_sources": {},
        "thresholds": {"boundary_distance_bp": BOUNDARY_BP, "af_common_threshold": 0.05, "min_sv_count_per_block": 3},
        "seeds": {"permutation_seed": 1, "umap_seed": 1},
        "paths": {
            "sv_calls": str(in_dir / "sv_calls.tsv"),
            "haploblocks": str(in_dir / "haploblocks.tsv"),
            "sample_metadata": str(in_dir / "sample_metadata.tsv"),
        },
    }
    with open(in_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(config, fh)

    out_dir = tmp_path_factory.mktemp("stage2_handcrafted_out")
    run(STAGE2_SCRIPT, ["--config", str(in_dir / "config.yaml"), "--out-dir", str(out_dir)])
    return pd.read_csv(out_dir / "sv_calls.tsv", sep="\t").set_index("sv_id")


def test_within_block(handcrafted_output):
    row = handcrafted_output.loc["sv_within"]
    assert row["position_class"] == "within_block"
    assert row["haploblock_id"] == "hbA"


def test_boundary_crossing_spans_two_blocks(handcrafted_output):
    row = handcrafted_output.loc["sv_span"]
    assert row["position_class"] == "boundary_crossing"
    assert set(row["haploblock_id"].split(",")) == {"hbA", "hbB"}


def test_boundary_crossing_near_shared_edge_lists_both_neighbors(handcrafted_output):
    row = handcrafted_output.loc["sv_near_shared_edge"]
    assert row["position_class"] == "boundary_crossing"
    assert set(row["haploblock_id"].split(",")) == {"hbA", "hbB"}


def test_boundary_crossing_near_isolated_edge_lists_one_block(handcrafted_output):
    row = handcrafted_output.loc["sv_near_isolated_edge"]
    assert row["position_class"] == "boundary_crossing"
    assert row["haploblock_id"] == "hbD"


def test_outside_block_far_from_blocks(handcrafted_output):
    row = handcrafted_output.loc["sv_far"]
    assert row["position_class"] == "outside_block"


def test_outside_block_no_haploblocks_on_chromosome(handcrafted_output):
    row = handcrafted_output.loc["sv_no_chrom"]
    assert row["position_class"] == "outside_block"


def test_overlapping_haploblocks_raise(tmp_path):
    haploblocks = pd.DataFrame(
        [
            {"haploblock_id": "hbA", "chrom": "chr1", "start": 1000, "end": 2000, "n_snps": 1, "n_clusters": 1, "hash_length": 1, "cluster_diff_score": 0.1},
            {"haploblock_id": "hbB", "chrom": "chr1", "start": 1500, "end": 2500, "n_snps": 1, "n_clusters": 1, "hash_length": 1, "cluster_diff_score": 0.1},
        ]
    )
    sv_calls = pd.DataFrame([{"sv_id": "sv1", "chrom": "chr1", "start": 1100, "end": 1200, "sv_type": "DEL", "imprecise": False, "length": 100}])
    sample_metadata = pd.DataFrame([{"sample_id": "SAMP0", "superpopulation": "EUR"}])
    haploblocks.to_csv(tmp_path / "haploblocks.tsv", sep="\t", index=False)
    sv_calls.to_csv(tmp_path / "sv_calls.tsv", sep="\t", index=False)
    sample_metadata.to_csv(tmp_path / "sample_metadata.tsv", sep="\t", index=False)
    config = {
        "thresholds": {"boundary_distance_bp": 100},
        "paths": {
            "sv_calls": str(tmp_path / "sv_calls.tsv"),
            "haploblocks": str(tmp_path / "haploblocks.tsv"),
            "sample_metadata": str(tmp_path / "sample_metadata.tsv"),
        },
    }
    with open(tmp_path / "config.yaml", "w") as fh:
        yaml.safe_dump(config, fh)

    result = subprocess.run(
        [sys.executable, str(STAGE2_SCRIPT), "--config", str(tmp_path / "config.yaml"), "--out-dir", str(tmp_path / "out")],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode != 0
    assert "overlap" in result.stderr


# --------------------------------------------------------------------------
# Integration: Stage 0 -> Stage 1 -> Stage 2 on synthetic (offline) data
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def stage2_from_pipeline(tmp_path_factory):
    stage0_dir = tmp_path_factory.mktemp("s0")
    run(
        STAGE0_SCRIPT,
        [
            "--skip-dbvar-download", "--skip-haploblock-download",
            "--out-dir", str(stage0_dir),
            "--n-svs", "200", "--n-haploblocks", "8", "--n-samples", "10",
            "--seed", "11",
        ],
    )
    stage1_dir = tmp_path_factory.mktemp("s1")
    run(STAGE1_SCRIPT, ["--config", str(stage0_dir / "config.yaml"), "--out-dir", str(stage1_dir)])
    stage2_dir = tmp_path_factory.mktemp("s2")
    run(STAGE2_SCRIPT, ["--config", str(stage1_dir / "config.yaml"), "--out-dir", str(stage2_dir)])
    return stage2_dir


def test_pipeline_wiring_produces_all_position_classes(stage2_from_pipeline):
    sv = pd.read_csv(stage2_from_pipeline / "sv_calls.tsv", sep="\t")
    assert len(sv) > 0
    assert "position_class" in sv.columns and "haploblock_id" in sv.columns
    assert set(sv["position_class"].unique()) <= {"within_block", "boundary_crossing", "outside_block"}
    # the synthetic generator places ~60/20/20% within/boundary/outside, so with
    # 200 SVs all three classes should show up
    assert set(sv["position_class"].unique()) == {"within_block", "boundary_crossing", "outside_block"}


def test_pipeline_stage2_config_points_at_its_own_output(stage2_from_pipeline):
    with open(stage2_from_pipeline / "config.yaml") as fh:
        config = yaml.safe_load(fh)
    for key in ("sv_calls", "haploblocks", "sample_metadata"):
        path = Path(config["paths"][key])
        assert path.parent == stage2_from_pipeline
        assert path.exists()
