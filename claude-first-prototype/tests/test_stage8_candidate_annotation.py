"""Tests for consequence-aware candidate annotation."""

import sys
from pathlib import Path

import pandas as pd
import yaml

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage8_candidate_annotation
import stage9_report


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
        "association_pattern": "cross_population_consistent_tag_candidate",
        "population_adjusted_r": 0.8,
        "q_value": 0.001, "informative_populations": 2, "directional_consistency": 1.0,
    }
    candidates = pd.DataFrame([
        {**common, "sv_id": "del", "start": 125, "end": 135, "sv_type": "DEL"},
        {**common, "sv_id": "del_exon", "start": 110, "end": 150, "sv_type": "DEL"},
        {**common, "sv_id": "ins", "start": 130, "end": 131, "sv_type": "INS"},
        {**common, "sv_id": "dup", "start": 90, "end": 210, "sv_type": "DUP"},
        {**common, "sv_id": "inv_break", "start": 150, "end": 280, "sv_type": "INV"},
        {**common, "sv_id": "inv_span", "start": 250, "end": 400, "sv_type": "INV"},
    ])

    serial_result = stage8_candidate_annotation.annotate_candidates(candidates, features)
    parallel_result = stage8_candidate_annotation.annotate_candidates(candidates, features, threads=2)
    pd.testing.assert_frame_equal(serial_result, parallel_result)
    result = serial_result.set_index("sv_id")
    assert result.loc["del", "consequence"] == "exonic_deletion"
    assert result.loc["del_exon", "consequence"] == "complete_exon_loss"
    assert result.loc["ins", "consequence"] == "exonic_insertion"
    assert result.loc["dup", "consequence"] == "complete_gene_duplication"
    assert result.loc["inv_break", "consequence"] == "inversion_breakpoint_disruption"
    assert result.loc["inv_span", "consequence"] == "gene_within_inversion"
    assert result.loc["inv_span", "overlap_basis"] == "span_only"
    assert set(result["call_quality"]) == {"pass_precise"}
    assert "candidate_score" not in result.columns


def test_representation_fields_are_carried_into_annotation() -> None:
    candidates = pd.DataFrame([{
        "sv_record_id": "record_1", "sv_id": "v1", "chrom": "chr1",
        "start": 10, "end": 11, "sv_type": "INS", "length": 1,
        "filter": "PASS", "imprecise": False, "haploblock_id": "block",
        "best_cluster_id": "C1", "association_pattern": "cluster_associated",
        "population_adjusted_r": 0.8, "q_value": 0.01,
        "informative_populations": 2, "directional_consistency": 1.0,
    }])
    representation = pd.DataFrame([{
        "sv_record_id": "record_1", "haploblock_id": "block",
        "representation_pattern": "multi_cluster_sv_candidate",
        "n_supported_carrier_clusters": 3,
        "top_cluster_carrier_evidence_share": 0.5,
        "top_standard_evidence_cluster_id": "C2",
        "top_standard_cluster_carrier_evidence_share": 0.8,
        "effective_standard_carrier_cluster_count": 1.5,
        "top_standard_cluster_carrier_rate": 0.9,
    }])
    result = stage8_candidate_annotation.annotate_candidates(
        candidates, pd.DataFrame(columns=["chrom", "feature", "start", "end", "gene_name"]),
        representation=representation,
    ).iloc[0]
    assert result["representation_pattern"] == "multi_cluster_sv_candidate"
    assert result["n_supported_carrier_clusters"] == 3
    assert result["top_standard_evidence_cluster_id"] == "C2"
    assert result["top_standard_cluster_carrier_evidence_share"] == 0.8


def test_stage8_triggers_report_after_writing_outputs(
    tmp_path: Path, monkeypatch
) -> None:
    candidates_path = tmp_path / "candidates.tsv"
    pd.DataFrame([{
        "sv_record_id": "record_1", "sv_id": "v1", "chrom": "chr1",
        "start": 10, "end": 20, "sv_type": "DEL", "length": 10,
        "filter": "PASS", "imprecise": False, "haploblock_id": "block",
        "best_cluster_id": "C1", "association_pattern": "cluster_associated",
        "population_adjusted_r": 0.8, "q_value": 0.01,
        "informative_populations": 2, "directional_consistency": 1.0,
    }]).to_csv(candidates_path, sep="\t", index=False)
    gtf_path = tmp_path / "genes.gtf"
    gtf_path.write_text(
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1"; gene_name "GENE1";\n'
    )
    config_path = tmp_path / "stage7.yaml"
    config_path.write_text(yaml.safe_dump({"paths": {
        "sv_cluster_summary": str(candidates_path), "gtf": str(gtf_path),
    }}))
    stage5_config = tmp_path / "stage5.yaml"
    stage5_config.write_text(yaml.safe_dump({"paths": {}}))
    captured = {}

    def capture_report(config_paths, out_dir, agent_mode, model) -> None:
        captured.update({
            "config_paths": config_paths, "out_dir": out_dir,
            "agent_mode": agent_mode, "model": model,
        })

    monkeypatch.setattr(stage9_report, "run_report", capture_report)
    stage8_out = tmp_path / "stage8"
    report_out = tmp_path / "stage9"
    stage8_candidate_annotation.main([
        "--config", str(config_path), "--out-dir", str(stage8_out),
        "--stage5-config", str(stage5_config), "--report-out-dir", str(report_out),
        "--report-agent", "off", "--report-model", "test-model",
    ])

    assert captured["config_paths"] == [stage5_config, stage8_out / "config.yaml"]
    assert captured["out_dir"] == report_out
    assert captured["agent_mode"] == "off"
    assert captured["model"] == "test-model"
    assert (stage8_out / "annotated_sv_candidates.tsv").exists()
