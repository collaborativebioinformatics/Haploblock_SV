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
            "sv_id": "portable", "chrom": "chr1", "start": 100, "end": 150,
            "sv_type": "DEL", "length": 50, "filter": "PASS", "imprecise": False,
        },
        {
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
        {"sv_id": "portable", "chrom": "chr1", "start": 100, "end": 150,
         "haploblock_id": "block1"},
        {"sv_id": "ancestry", "chrom": "chr1", "start": 200, "end": 250,
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

    portable = summary[summary["sv_id"] == "portable"].iloc[0]
    assert portable["association_pattern"] == "portable_cluster_tag"
    ancestry = summary[summary["sv_id"] == "ancestry"].iloc[0]
    assert ancestry["association_pattern"] == "no_detected_cluster_signal"


def test_variant_coordinates_disambiguate_reused_vcf_ids() -> None:
    samples = [f"S{i}" for i in range(8)]
    metadata = pd.DataFrame({"sample_id": samples, "population": ["P"] * 8})
    rows = []
    for start, carriers in ((100, {"S0", "S1", "S2", "S3"}), (200, set())):
        row = {
            "sv_id": "reused", "chrom": "chr1", "start": start, "end": start + 10,
            "sv_type": "INS", "length": 10, "filter": "PASS", "imprecise": False,
        }
        row.update({sample: "0/1" if sample in carriers else "0/0" for sample in samples})
        rows.append(row)
    blocks = pd.DataFrame([
        {"sv_id": "reused", "chrom": "chr1", "start": start, "end": start + 10,
         "haploblock_id": "block"}
        for start in (100, 200)
    ])
    memberships = pd.DataFrame([
        {"haploblock_id": "block", "sample_id": sample, "haplotype": haplotype,
         "cluster_id": "tag" if int(sample[1:]) < 4 else "other"}
        for sample in samples for haplotype in (0, 1)
    ])
    result = stage6_cluster_association.association_table(
        pd.DataFrame(rows), blocks, memberships, metadata,
        permutations=19, seed=2, min_cluster_haplotypes=4, min_population_samples=4,
    )
    assert set(result["start"]) == {100}

