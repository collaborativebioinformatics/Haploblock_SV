"""Tests for population-conditioned SV-cluster association."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage6_cluster_association


def test_vectorized_population_adjustment_matches_scalar_correlations() -> None:
    populations = np.array(["A", "A", "A", "B", "B", "B"])
    cluster_matrix = np.array([
        [0, 2], [1, 1], [2, 0], [0, 1], [1, 0], [2, 2],
    ], dtype=float)
    y = np.array([
        [0, 1, 2, 0, 1, 2],
        [0, np.nan, 1, 2, 1, 0],
    ], dtype=float)
    called = np.isfinite(y)
    observed, _ = stage6_cluster_association.population_adjusted_correlations(
        y, called, cluster_matrix, populations
    )

    for variant_index in range(len(y)):
        use = called[variant_index]
        for cluster_index in range(cluster_matrix.shape[1]):
            expected = stage6_cluster_association.correlation(
                cluster_matrix[use, cluster_index],
                y[variant_index, use],
                populations[use],
            )
            assert np.isclose(observed[variant_index, cluster_index], expected)


def test_staged_permutations_only_advance_selected_variants() -> None:
    populations = np.array(["P"] * 6)
    cluster_matrix = np.array([
        [0, 1], [0, 2], [1, 1], [1, 0], [2, 0], [2, 2],
    ], dtype=float)
    y = np.array([
        [0, 0, 1, 1, 2, 2],
        [0, 1, 0, 2, 1, 2],
    ], dtype=float)
    called = np.isfinite(y)
    observed, cluster_variance = (
        stage6_cluster_association.population_adjusted_correlations(
            y, called, cluster_matrix, populations
        )
    )
    eligible = cluster_variance > 0

    _, triage_only = stage6_cluster_association.staged_permutation_p_values(
        y, called, cluster_matrix, populations, observed, eligible,
        permutations=19, refinement_permutations=None,
        refinement_p_threshold=0.01, seed=8,
        initial_permutations=5, initial_p_threshold=0.0,
    )
    assert set(triage_only) == {5}

    _, screened = stage6_cluster_association.staged_permutation_p_values(
        y, called, cluster_matrix, populations, observed, eligible,
        permutations=19, refinement_permutations=None,
        refinement_p_threshold=0.01, seed=8,
        initial_permutations=5, initial_p_threshold=1.0,
    )
    assert set(screened) == {19}

    _, effect_filtered = stage6_cluster_association.staged_permutation_p_values(
        y, called, cluster_matrix, populations, observed, eligible,
        permutations=19, refinement_permutations=None,
        refinement_p_threshold=0.01, seed=8,
        initial_permutations=5, initial_p_threshold=1.0,
        min_refinement_abs_r=2.0,
    )
    assert set(effect_filtered) == {5}


def test_within_population_test_separates_portable_from_ancestry_signal() -> None:
    samples = [f"A{i}" for i in range(8)] + [f"B{i}" for i in range(8)]
    metadata = pd.DataFrame({
        "sample_id": samples,
        "population": ["A"] * 8 + ["B"] * 8,
    })
    variant_rows = [
        {
            "sv_record_id": "record_portable",
            "sv_id": "portable", "chrom": "chr1", "start": 100, "end": 150,
            "sv_type": "DEL", "length": 50, "filter": "PASS", "imprecise": False,
        },
        {
            "sv_record_id": "record_ancestry",
            "sv_id": "ancestry", "chrom": "chr1", "start": 200, "end": 250,
            "sv_type": "DEL", "length": 50, "filter": "PASS", "imprecise": False,
        },
    ]
    for sample in samples:
        index = int(sample[1:])
        variant_rows[0][sample] = "0/1" if index < 4 else "0/0"
        variant_rows[1][sample] = "0/1" if sample.startswith("A") else "0/0"
    sv = pd.DataFrame(variant_rows)
    sv_blocks = pd.DataFrame([
        {"sv_record_id": "record_portable", "sv_id": "portable", "chrom": "chr1", "start": 100, "end": 150,
         "haploblock_id": "block1"},
        {"sv_record_id": "record_ancestry", "sv_id": "ancestry", "chrom": "chr1", "start": 200, "end": 250,
         "haploblock_id": "block2"},
    ])

    memberships = []
    for sample in samples:
        index = int(sample[1:])
        local_cluster = "tag" if index < 4 else "other"
        ancestry_cluster = "A_cluster" if sample.startswith("A") else "B_cluster"
        for haplotype in (0, 1):
            memberships.extend([
                {"haploblock_id": "block1", "sample_id": sample,
                 "haplotype": haplotype, "cluster_id": local_cluster},
                {"haploblock_id": "block2", "sample_id": sample,
                 "haplotype": haplotype, "cluster_id": ancestry_cluster},
            ])

    result = stage6_cluster_association.association_table(
        sv, sv_blocks, pd.DataFrame(memberships), metadata,
        permutations=199, seed=4, min_cluster_haplotypes=4, min_population_samples=4,
    )
    parallel_result = stage6_cluster_association.association_table(
        sv, sv_blocks, pd.DataFrame(memberships), metadata,
        permutations=199, seed=4, min_cluster_haplotypes=4, min_population_samples=4,
        threads=2,
    )
    pd.testing.assert_frame_equal(result, parallel_result)
    summary = stage6_cluster_association.summarize_associations(
        result, q_threshold=0.05, min_abs_r=0.3
    )
    assert result.groupby(["sv_record_id", "haploblock_id"])["q_value"].nunique().max() == 1

    portable = summary[summary["sv_id"] == "portable"].iloc[0]
    assert portable["association_pattern"] == "cross_population_consistent_tag_candidate"
    assert portable["association_direction"] == "carrier_enriched"
    tag = result[(result["sv_id"] == "portable") & (result["cluster_id"] == "tag")].iloc[0]
    assert tag["n_samples_with_cluster"] == 8
    assert tag["n_sv_carriers_with_cluster"] == 8
    assert tag["n_sv_noncarriers_with_cluster"] == 0
    ancestry = summary[summary["sv_id"] == "ancestry"].iloc[0]
    assert ancestry["association_pattern"] == "no_detected_cluster_signal"


def test_summary_keeps_strong_exclusion_over_weak_enrichment() -> None:
    common = {
        "sv_record_id": "record", "sv_id": "v", "chrom": "chr1", "start": 1,
        "end": 2, "sv_type": "INS", "length": 1, "filter": "PASS",
        "imprecise": False, "haploblock_id": "block", "n_called": 20,
        "n_cluster_carriers": 10, "n_samples_with_cluster": 10,
        "n_sv_carriers_with_cluster": 5, "n_sv_noncarriers_with_cluster": 5,
        "cluster_haplotype_count": 10, "carrier_rate_with_cluster": 0.5,
        "carrier_rate_without_cluster": 0.5, "informative_populations": 2,
        "directional_consistency": 1.0, "q_value": 0.01,
    }
    associations = pd.DataFrame([
        {**common, "cluster_id": "weak_positive", "carrier_rate_difference": 0.05,
         "association_direction": "carrier_enriched", "population_adjusted_r": 0.1,
         "p_value": 0.2},
        {**common, "cluster_id": "strong_negative", "carrier_rate_difference": -0.8,
         "association_direction": "carrier_depleted", "population_adjusted_r": -0.9,
         "p_value": 0.01},
    ])
    summary = stage6_cluster_association.summarize_associations(
        associations, q_threshold=0.05, min_abs_r=0.3
    ).iloc[0]
    assert summary["best_cluster_id"] == "strong_negative"
    assert summary["best_enriched_cluster_id"] == "weak_positive"
    assert summary["association_pattern"] == "cluster_exclusion_signal"


def test_cluster_present_in_every_sample_has_no_exclusion_comparator() -> None:
    samples = [f"S{i}" for i in range(8)]
    metadata = pd.DataFrame({"sample_id": samples, "population": ["P"] * 8})
    variant = {
        "sv_record_id": "record", "sv_id": "v", "chrom": "chr1", "start": 1,
        "end": 2, "sv_type": "INS", "length": 1, "filter": "PASS", "imprecise": False,
    }
    variant.update({sample: "0/1" if index < 4 else "0/0" for index, sample in enumerate(samples)})
    blocks = pd.DataFrame([{
        "sv_record_id": "record", "sv_id": "v", "chrom": "chr1", "start": 1,
        "end": 2, "haploblock_id": "block",
    }])
    memberships = pd.DataFrame([
        {
            "haploblock_id": "block", "sample_id": sample, "haplotype": haplotype,
            "cluster_id": "tag" if index < 4 or haplotype == 0 else "other",
        }
        for index, sample in enumerate(samples) for haplotype in (0, 1)
    ])

    result = stage6_cluster_association.association_table(
        pd.DataFrame([variant]), blocks, memberships, metadata,
        permutations=19, seed=3, min_cluster_haplotypes=4, min_population_samples=4,
    )
    tag = result[result["cluster_id"] == "tag"].iloc[0]
    assert pd.isna(tag["carrier_rate_without_cluster"])
    assert pd.isna(tag["carrier_rate_difference"])
    assert tag["association_direction"] == "unavailable_comparator"

    summary = stage6_cluster_association.summarize_associations(
        result, q_threshold=1.1, min_abs_r=0.0
    )
    assert summary.iloc[0]["association_pattern"] != "cluster_exclusion_signal"


def test_directional_consistency_uses_mean_dosage_not_carrier_rate() -> None:
    samples = [f"S{i}" for i in range(6)]
    metadata = pd.DataFrame({"sample_id": samples, "population": ["P"] * 6})
    variant = {
        "sv_record_id": "record", "sv_id": "v", "chrom": "chr1", "start": 1,
        "end": 2, "sv_type": "INS", "length": 1, "filter": "PASS", "imprecise": False,
    }
    dosages = ["0/1", "0/1", "0/1", "1/1", "1/1", "0/0"]
    variant.update(dict(zip(samples, dosages)))
    blocks = pd.DataFrame([{
        "sv_record_id": "record", "sv_id": "v", "chrom": "chr1", "start": 1,
        "end": 2, "haploblock_id": "block",
    }])
    memberships = pd.DataFrame([
        {
            "haploblock_id": "block", "sample_id": sample, "haplotype": haplotype,
            "cluster_id": "tag" if index < 3 else "other",
        }
        for index, sample in enumerate(samples) for haplotype in (0, 1)
    ])

    result = stage6_cluster_association.association_table(
        pd.DataFrame([variant]), blocks, memberships, metadata,
        permutations=19, seed=5, min_cluster_haplotypes=4, min_population_samples=4,
    )
    tag = result[result["cluster_id"] == "tag"].iloc[0]
    assert tag["carrier_rate_difference"] > 0
    assert tag["population_adjusted_r"] < 0
    assert tag["directional_consistency"] == 1.0


def test_record_ids_disambiguate_reused_vcf_ids_and_coordinates() -> None:
    samples = [f"S{i}" for i in range(8)]
    metadata = pd.DataFrame({"sample_id": samples, "population": ["P"] * 8})
    rows = []
    for record_id, carriers in (("record_1", {"S0", "S1", "S2", "S3"}), ("record_2", set())):
        row = {
            "sv_record_id": record_id,
            "sv_id": "reused", "chrom": "chr1", "start": 100, "end": 110,
            "sv_type": "INS", "length": 10, "filter": "PASS", "imprecise": False,
        }
        row.update({sample: "0/1" if sample in carriers else "0/0" for sample in samples})
        rows.append(row)
    blocks = pd.DataFrame([
        {"sv_record_id": record_id, "sv_id": "reused", "chrom": "chr1", "start": 100, "end": 110,
             "haploblock_id": "block"}
        for record_id in ("record_1", "record_2")
    ])
    memberships = pd.DataFrame([
        {"haploblock_id": "block", "sample_id": sample, "haplotype": haplotype,
         "cluster_id": "tag" if int(sample[1:]) < 4 else "other"}
        for sample in samples for haplotype in (0, 1)
    ])
    result = stage6_cluster_association.association_table(
        pd.DataFrame(rows), blocks, memberships, metadata,
        permutations=19, seed=2, min_cluster_haplotypes=4, min_population_samples=4,
        refinement_permutations=29, refinement_p_threshold=1.0,
    )
    assert set(result["sv_record_id"]) == {"record_1"}
    assert set(result["permutations_used"]) == {29}
