"""Tests for haploblock information gained from SVs."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage7_information_gain


def test_dosage_fast_path_falls_back_for_multiallelic_genotypes() -> None:
    genotypes = pd.DataFrame({
        "S1": ["0|1", "./."],
        "S2": ["0|2", "2|2"],
    })
    result = stage7_information_gain.dosage_matrix_from_genotypes(
        genotypes, ["S1", "S2"]
    )
    np.testing.assert_equal(result, np.array([[1.0, 1.0], [np.nan, 2.0]]))


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

    dosages = stage7_information_gain.dosage_table(sv, samples)
    reused = stage7_information_gain.information_table(
        sv, blocks, memberships, samples, min_diplotype_samples=4, dosages=dosages
    )
    pd.testing.assert_frame_equal(
        result.reset_index().sort_index(axis=1), reused.sort_index(axis=1)
    )

    without_secondary_metric = stage7_information_gain.information_table(
        sv, blocks, memberships, samples, min_diplotype_samples=4,
        dosages=dosages, include_information_gain=False,
    ).set_index("sv_id")
    assert without_secondary_metric["normalized_information_gain"].isna().all()
    assert without_secondary_metric.loc["mixed", "mixed_diplotype_fraction"] == 1.0

    assignments = pd.DataFrame([{
        **{column: tagged[column] for column in stage7_information_gain.METADATA_COLUMNS},
        "haploblock_id": "block", "cluster_id": "C1",
        "expected_alt_haplotypes": 8.0, "evidence_tier": "standard",
    }])
    carrier_clusters = stage7_information_gain.carrier_cluster_summary(assignments)
    purity = stage7_information_gain.cluster_purity_table(
        sv, blocks, memberships, samples, min_cluster_samples=4,
        min_carriers=3, min_noncarriers=3,
    )
    representation = stage7_information_gain.representation_summary(
        carrier_clusters, purity, result.reset_index(), purity_threshold=0.9
    ).set_index("sv_id")
    assert representation.loc["tagged", "n_supported_carrier_clusters"] == 1
    assert representation.loc["tagged", "top_cluster_carrier_evidence_share"] == 1.0
    assert representation.loc["tagged", "representation_pattern"] == "hash_tag_candidate"
    assert representation.loc["mixed", "n_mixed_diplotypes_meeting_count_threshold"] == 2
    assert representation.loc["mixed", "representation_pattern"] == "hash_subdivision_candidate"


def test_carrier_cluster_summary_detects_multi_cluster_sv() -> None:
    base = {
        "sv_record_id": "record_multi", "sv_id": "multi", "chrom": "chr1",
        "start": 100, "end": 110, "sv_type": "INS", "length": 10,
        "filter": "PASS", "imprecise": False, "haploblock_id": "block",
        "evidence_tier": "standard",
    }
    assignments = pd.DataFrame([
        {**base, "cluster_id": "C1", "expected_alt_haplotypes": 6.0},
        {**base, "cluster_id": "C2", "expected_alt_haplotypes": 2.0},
    ])
    summary = stage7_information_gain.carrier_cluster_summary(assignments).iloc[0]
    assert summary["n_supported_carrier_clusters"] == 2
    assert summary["top_cluster_carrier_evidence_share"] == 0.75
    assert summary["effective_carrier_cluster_count"] == 1.6

    mixed_evidence = assignments.copy()
    mixed_evidence.loc[mixed_evidence["cluster_id"] == "C1", "evidence_tier"] = "low"
    mixed_summary = stage7_information_gain.carrier_cluster_summary(mixed_evidence).iloc[0]
    assert mixed_summary["top_supported_cluster_id"] == "C1"
    assert mixed_summary["top_standard_evidence_cluster_id"] == "C2"

    representation = stage7_information_gain.representation_summary(
        stage7_information_gain.carrier_cluster_summary(
            assignments.assign(evidence_tier="low")
        ),
        pd.DataFrame(), pd.DataFrame(), 0.9,
    ).iloc[0]
    assert representation["n_supported_carrier_clusters"] == 2
    assert representation["n_standard_evidence_carrier_clusters"] == 0
    assert representation["representation_pattern"] == "insufficient_or_partial_evidence"


def test_population_specific_sv_on_shared_cluster_is_flagged() -> None:
    metadata = {
        "sv_record_id": "record", "sv_id": "v", "chrom": "chr1", "start": 1,
        "end": 2, "sv_type": "INS", "length": 1, "filter": "PASS", "imprecise": False,
    }
    carrier_clusters = pd.DataFrame([{
        **metadata, "haploblock_id": "block", "n_supported_carrier_clusters": 1,
        "n_standard_evidence_carrier_clusters": 1, "top_supported_cluster_id": "C1",
        "top_cluster_carrier_evidence_share": 1.0, "effective_carrier_cluster_count": 1.0,
        "top_standard_evidence_cluster_id": "C1",
        "top_standard_cluster_carrier_evidence_share": 1.0,
        "effective_standard_carrier_cluster_count": 1.0,
    }])
    purity = pd.DataFrame([{
        **metadata, "haploblock_id": "block", "cluster_id": "C1",
        "n_called_cluster_samples": 10, "n_sv_carriers": 10, "n_sv_noncarriers": 0,
        "carrier_rate_in_cluster": 1.0, "cluster_purity": 1.0, "mixed_balance": 0.0,
        "meets_mixed_count_threshold": False,
    }])
    cluster_populations = pd.DataFrame([{
        "haploblock_id": "block", "cluster_id": "C1",
        "cluster_population_count": 3, "cluster_populations": "A;B;C",
    }])
    classifications = pd.DataFrame([{
        "sv_record_id": "record", "sv_class": "population_specific",
        "specific_to_population": "A",
    }])
    result = stage7_information_gain.representation_summary(
        carrier_clusters, purity, pd.DataFrame(), 0.9,
        cluster_populations, classifications, 3,
    ).iloc[0]
    assert result["population_context_pattern"] == (
        "population_enriched_on_shared_cluster_candidate"
    )


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
    assert np.isclose(variance["explained_variance_ratio"].sum(), 1.0)


def test_stage7_speed_flags_are_optional() -> None:
    defaults = stage7_information_gain.parse_args(["--config", "config.yaml"])
    assert defaults.threads == 8
    assert not defaults.skip_information_gain
    assert not defaults.skip_pca

    accelerated = stage7_information_gain.parse_args([
        "--config", "config.yaml", "--threads", "3",
        "--skip-information-gain", "--skip-pca",
    ])
    assert accelerated.threads == 3
    assert accelerated.skip_information_gain
    assert accelerated.skip_pca
