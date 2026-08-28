"""Tests for population-conditioned SV-cluster association."""

import sys
from pathlib import Path

import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage6_cluster_association


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
