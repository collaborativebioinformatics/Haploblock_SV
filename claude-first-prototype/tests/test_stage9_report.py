"""Tests for Stage 9 deterministic summaries and agent handoff."""

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage9_report


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_report_writes_traceable_facts_and_figures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    enrichment_path = tmp_path / "enrichment.tsv"
    associations_path = tmp_path / "associations.tsv"
    representation_path = tmp_path / "representation.tsv"
    annotations_path = tmp_path / "annotations.tsv"
    write_table(enrichment_path, [
        {"haploblock_id": "b1", "sv_type": "DEL", "observed_count": 8,
         "expected_count": 2.0, "p_value": 0.001, "q_value": 0.01, "flagged": True},
        {"haploblock_id": "b2", "sv_type": "INS", "observed_count": 1,
         "expected_count": 1.0, "p_value": 1.0, "q_value": 1.0, "flagged": False},
    ])
    write_table(associations_path, [
        {"sv_record_id": "r1", "sv_id": "v1", "chrom": "chr1", "start": 10,
         "end": 20, "sv_type": "DEL", "haploblock_id": "b1",
         "best_cluster_id": "c1", "population_adjusted_r": 0.8,
         "p_value": 0.001, "q_value": 0.02, "carrier_rate_with_cluster": 0.9,
         "carrier_rate_without_cluster": 0.1,
         "permutations_used": 10000,
         "association_pattern": "cross_population_consistent_tag_candidate"},
        {"sv_record_id": "r2", "sv_id": "v2", "chrom": "chr1", "start": 30,
         "end": 31, "sv_type": "INS", "haploblock_id": "b2",
         "best_cluster_id": "c2", "population_adjusted_r": 0.1,
         "p_value": 0.5, "q_value": 0.8, "carrier_rate_with_cluster": 0.2,
         "carrier_rate_without_cluster": 0.1,
         "permutations_used": 200,
         "association_pattern": "no_detected_cluster_signal"},
    ])
    write_table(representation_path, [
        {"sv_record_id": "r1", "haploblock_id": "b1",
         "representation_pattern": "hash_tag_candidate",
         "population_context_pattern": "no_population_restriction_pattern"},
    ])
    write_table(annotations_path, [
        {"sv_record_id": "r1", "sv_id": "v1", "chrom": "chr1", "start": 10,
         "end": 20, "sv_type": "DEL", "haploblock_id": "b1", "genes": "GENE1",
         "consequence": "complete_exon_loss", "overlap_basis": "exon",
         "filter": "PASS", "imprecise": False, "population_adjusted_r": 0.8,
         "q_value": 0.02,
         "association_pattern": "cross_population_consistent_tag_candidate",
         "representation_pattern": "hash_tag_candidate", "sv_class": "common",
         "specific_to_population": None},
        {"sv_record_id": "r2", "sv_id": "v2", "chrom": "chr1", "start": 30,
         "end": 31, "sv_type": "INS", "haploblock_id": "b2", "genes": "GENE2",
         "consequence": "complete_exon_loss", "overlap_basis": "exon",
         "filter": "PASS", "imprecise": False, "population_adjusted_r": 0.9,
         "q_value": 0.02,
         "association_pattern": "cross_population_consistent_tag_candidate",
         "representation_pattern": "not_evaluated", "sv_class": "common",
         "specific_to_population": None},
    ])
    stage5_config = tmp_path / "stage5.yaml"
    stage8_config = tmp_path / "stage8.yaml"
    stage5_config.write_text(yaml.safe_dump({"paths": {
        "sv_type_enrichment": str(enrichment_path),
    }}))
    stage8_config.write_text(yaml.safe_dump({"paths": {
        "sv_cluster_summary": str(associations_path),
        "sv_hash_representation": str(representation_path),
        "annotated_sv_candidates": str(annotations_path),
    }}))

    out_dir = tmp_path / "report"
    monkeypatch.setenv("OPENAI_API_KEY", "test-key")

    def fail_agent(*args: object, **kwargs: object) -> tuple[str, str | None]:
        raise RuntimeError("temporary API failure")

    monkeypatch.setattr(stage9_report, "request_agent_interpretation", fail_agent)
    stage9_report.main([
        "--config", str(stage5_config), "--config", str(stage8_config),
        "--out-dir", str(out_dir), "--agent", "auto",
    ])

    facts = json.loads((out_dir / "report_facts.json").read_text())
    assert facts["stage5_sv_type_enrichment"]["n_significant_q_lt_0_05"] == 1
    assert facts["stage5_sv_type_enrichment"]["top_cells"][0]["flagged"] is True
    assert facts["stage6_cluster_association"]["n_significant_q_lt_0_05"] == 1
    assert facts["stage6_cluster_association"]["permutations_used_counts"] == {
        "10000": 1, "200": 1,
    }
    assert facts["stage7_hash_representation_and_qc"][
        "representation_pattern_counts"
    ] == {"hash_tag_candidate": 1}
    assert facts["stage8_candidate_annotation"]["top_candidates"][0]["genes"] == "GENE1"
    assert (out_dir / "report.md").exists()
    assert (out_dir / "report.html").exists()
    assert (out_dir / "figures/stage5_type_enrichment.png").exists()
    assert (out_dir / "figures/stage6_association_patterns.png").exists()
    assert (out_dir / "figures/stage7_representation_patterns.png").exists()
    assert (out_dir / "figures/stage8_consequences.png").exists()
    metadata = json.loads((out_dir / "agent_metadata.json").read_text())
    assert metadata["agent_used"] is False
    assert "temporary API failure" in metadata["error"]

    monkeypatch.delenv("OPENAI_API_KEY")
    monkeypatch.setattr(stage9_report, "create_figures", lambda *args: [])
    required_dir = tmp_path / "required_report"
    with pytest.raises(RuntimeError, match="OPENAI_API_KEY"):
        stage9_report.run_report(
            [stage5_config, stage8_config], required_dir, agent_mode="required"
        )
    assert (required_dir / "report.md").exists()
    assert (required_dir / "report.html").exists()
    assert (required_dir / "agent_metadata.json").exists()


def test_agent_receives_only_serialized_facts() -> None:
    captured = {}

    class FakeResponses:
        def create(self, **kwargs: object) -> SimpleNamespace:
            captured.update(kwargs)
            return SimpleNamespace(output_text="## Main findings\nA careful result.", id="resp_1")

    client = SimpleNamespace(responses=FakeResponses())
    text, response_id = stage9_report.request_agent_interpretation(
        {"stage": {"count": 2}}, "test-model", client
    )

    assert text.endswith("A careful result.")
    assert response_id == "resp_1"
    assert captured["model"] == "test-model"
    assert json.loads(captured["input"])["stage"]["count"] == 2
    assert "Every number must come directly" in captured["instructions"]
