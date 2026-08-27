"""Tests for Stage 4 population allele-frequency classification."""

import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage4_classify_af


def test_stage4_reads_stage1_contract_and_classifies_each_sv_once(tmp_path: Path) -> None:
    metadata_path = tmp_path / "sample_metadata.tsv"
    pd.DataFrame(
        {
            "sample_id": ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"],
            "population": ["popA"] * 4 + ["popB"] * 4,
        }
    ).to_csv(metadata_path, sep="\t", index=False)

    metadata_columns = {
        "sv_record_id": "",
        "chrom": "chr1",
        "start": 100,
        "end": 200,
        "sv_type": "DEL",
        "length": 100,
        "filter": "PASS",
        "imprecise": False,
    }
    sample_ids = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"]

    def variant(sv_id: str, genotypes: list[str]) -> dict:
        return {
            **metadata_columns, "sv_id": sv_id, "sv_record_id": f"record_{sv_id}",
            **dict(zip(sample_ids, genotypes)),
        }

    genotype_path = tmp_path / "sv_genotypes.chr1.tsv"
    pd.DataFrame(
        [
            variant("common", ["0/1", "0/0", "0/0", "0/0", "0/1", "0/0", "0/0", "0/0"]),
            variant("private", ["1/1", "0/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]),
            variant("rare", ["0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/."]),
        ]
    ).to_csv(genotype_path, sep="\t", index=False)

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "genome_build": "GRCh38",
                "thresholds": {"association_probability": 0.75},
                "paths": {
                    "sample_metadata": str(metadata_path),
                    "sv_genotypes": {"chr1": str(genotype_path)},
                },
            },
            sort_keys=False,
        )
    )

    out_dir = tmp_path / "stage4_output"
    stage4_classify_af.main(["--config", str(config_path), "--out-dir", str(out_dir)])

    by_population = pd.read_csv(out_dir / "sv_af_classification.tsv", sep="\t")
    classifications = pd.read_csv(out_dir / "sv_classification.tsv", sep="\t")
    assert len(by_population) == 3 * 2
    assert len(classifications) == 3
    assert "cluster_id" not in by_population.columns

    classes = classifications.set_index("sv_id")
    assert classes.loc["common", "sv_class"] == "common"
    assert classes.loc["private", "sv_class"] == "population_specific"
    assert classes.loc["private", "specific_to_population"] == "popA"
    assert classes.loc["rare", "sv_class"] == "other"
    assert classes.loc["rare", "other_reason"] == "absent_or_rare"

    af = by_population.set_index(["sv_id", "population"])["af"]
    assert af[("common", "popA")] == pytest.approx(1 / 8)
    assert af[("common", "popB")] == pytest.approx(1 / 8)
    assert af[("private", "popA")] == pytest.approx(3 / 8)
    assert af[("private", "popB")] == pytest.approx(0)

    rare_b = by_population[
        (by_population["sv_id"] == "rare") & (by_population["population"] == "popB")
    ].iloc[0]
    assert rare_b["n_called"] == 3
    assert rare_b["called_alleles"] == 6

    output_config = yaml.safe_load((out_dir / "config.yaml").read_text())
    assert output_config["paths"]["sv_classification"] == str(
        (out_dir / "sv_classification.tsv").resolve()
    )


def test_genotype_counts_treat_multidigit_alleles_as_one_allele() -> None:
    alternate, called = stage4_classify_af.genotype_counts(pd.Series(["12/0", "10/10", "./."]))

    assert alternate.tolist() == [1, 2, 0]
    assert called.tolist() == [2, 2, 0]
