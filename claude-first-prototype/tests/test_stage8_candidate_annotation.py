"""Tests for consequence-aware candidate annotation."""

import sys
from pathlib import Path

import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage8_candidate_annotation


def test_sv_types_receive_distinct_gene_consequences() -> None:
    features = pd.DataFrame([
        {"chrom": "chr1", "feature": "gene", "start": 100, "end": 200,
         "gene_id": "g1", "gene_name": "GENE1", "gene_biotype": "protein_coding"},
        {"chrom": "chr1", "feature": "exon", "start": 120, "end": 140,
         "gene_id": "g1", "gene_name": "GENE1", "gene_biotype": "protein_coding"},
        {"chrom": "chr1", "feature": "gene", "start": 300, "end": 350,
         "gene_id": "g2", "gene_name": "GENE2", "gene_biotype": "protein_coding"},
    ])
    common = {
        "chrom": "chr1", "length": 20, "filter": "PASS", "imprecise": False,
        "haploblock_id": "b", "best_cluster_id": "c",
        "association_pattern": "portable_cluster_tag", "population_adjusted_r": 0.8,
        "q_value": 0.001, "informative_populations": 2, "directional_consistency": 1.0,
    }
    candidates = pd.DataFrame([
        {**common, "sv_id": "del", "start": 125, "end": 135, "sv_type": "DEL"},
        {**common, "sv_id": "ins", "start": 130, "end": 131, "sv_type": "INS"},
        {**common, "sv_id": "dup", "start": 90, "end": 210, "sv_type": "DUP"},
        {**common, "sv_id": "inv_break", "start": 150, "end": 280, "sv_type": "INV"},
        {**common, "sv_id": "inv_span", "start": 250, "end": 400, "sv_type": "INV"},
    ])

    result = stage8_candidate_annotation.annotate_candidates(candidates, features).set_index("sv_id")
    assert result.loc["del", "consequence"] == "exon_loss"
    assert result.loc["ins", "consequence"] == "exonic_insertion"
    assert result.loc["dup", "consequence"] == "complete_gene_duplication"
    assert result.loc["inv_break", "consequence"] == "inversion_breakpoint_disruption"
    assert result.loc["inv_span", "consequence"] == "gene_within_inversion"
    assert result.loc["inv_span", "overlap_basis"] == "span_only"

