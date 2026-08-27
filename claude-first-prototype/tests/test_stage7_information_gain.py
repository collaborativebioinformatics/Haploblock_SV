"""Tests for haploblock information gained from SVs."""

import sys
from pathlib import Path

import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage7_information_gain


def test_information_gain_distinguishes_tagged_and_cluster_subdividing_svs() -> None:
    samples = [f"S{i}" for i in range(16)]
    base = {
        "chrom": "chr1", "sv_type": "DEL", "length": 10,
        "filter": "PASS", "imprecise": False,
    }
    tagged = {**base, "sv_record_id": "record_tagged", "sv_id": "tagged", "start": 100, "end": 110}
    mixed = {**base, "sv_record_id": "record_mixed", "sv_id": "mixed", "start": 200, "end": 210}
    for index, sample in enumerate(samples):
        tagged[sample] = "0/1" if index < 8 else "0/0"
        mixed[sample] = "0/1" if index % 2 == 0 else "0/0"
    sv = pd.DataFrame([tagged, mixed])
    blocks = pd.DataFrame([
            {"sv_record_id": row["sv_record_id"], "sv_id": row["sv_id"], "chrom": "chr1", "start": row["start"],
         "end": row["end"], "haploblock_id": "block"}
        for row in (tagged, mixed)
    ])
    memberships = pd.DataFrame([
        {"haploblock_id": "block", "sample_id": sample, "haplotype": haplotype,
         "cluster_id": "C1" if index < 8 else "C2"}
        for index, sample in enumerate(samples) for haplotype in (0, 1)
    ])

    result = stage7_information_gain.information_table(
        sv, blocks, memberships, samples, min_diplotype_samples=4
    ).set_index("sv_id")
    assert result.loc["tagged", "normalized_information_gain"] == 1.0
    assert result.loc["tagged", "mixed_diplotype_fraction"] == 0.0
    assert result.loc["mixed", "normalized_information_gain"] == 0.0
    assert result.loc["mixed", "mixed_diplotype_fraction"] == 1.0


def test_pca_writes_reusable_deterministic_tables() -> None:
    samples = ["S1", "S2", "S3", "S4"]
    metadata = pd.DataFrame({
        "sample_id": samples,
        "population": ["A", "A", "B", "B"],
        "superpopulation": ["X", "X", "Y", "Y"],
    })
    sv = pd.DataFrame([
        {"sv_id": "v1", "chrom": "chr1", "start": 1, "end": 2, "sv_type": "INS",
         "length": 1, "filter": "PASS", "imprecise": False,
         "S1": "0/0", "S2": "0/0", "S3": "1/1", "S4": "1/1"},
        {"sv_id": "v2", "chrom": "chr1", "start": 3, "end": 4, "sv_type": "INS",
         "length": 1, "filter": "PASS", "imprecise": False,
         "S1": "0/0", "S2": "0/1", "S3": "0/1", "S4": "1/1"},
    ])
    first, variance = stage7_information_gain.pca_tables(sv, metadata, 1.0, 0.01, 2)
    second, _ = stage7_information_gain.pca_tables(sv, metadata, 1.0, 0.01, 2)
    pd.testing.assert_frame_equal(first, second)
    assert list(variance["n_variants"]) == [2, 2]
