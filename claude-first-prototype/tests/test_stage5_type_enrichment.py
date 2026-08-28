"""Tests for Stage 5's reconciled SV-block input contract."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage5_type_enrichment


def test_stage5_uses_unique_sv_block_pairs_and_flags_spike(tmp_path: Path) -> None:
    haploblocks = pd.DataFrame(
        [
            {"haploblock_id": "small", "chrom": "chr1", "start": 0, "end": 1_000},
            {"haploblock_id": "medium", "chrom": "chr1", "start": 1_000, "end": 21_000},
            {"haploblock_id": "large", "chrom": "chr1", "start": 21_000, "end": 221_000},
        ]
    )
    rows = []

    def add_variants(block: str, sv_type: str, count: int) -> None:
        for _ in range(count):
            rows.append(
                    {
                        "sv_record_id": f"record{len(rows)}",
                        "sv_id": f"sv{len(rows)}",
                    "haploblock_id": block,
                    "chrom": "chr1",
                    "sv_type": sv_type,
                    "start": len(rows) * 100,
                    "end": len(rows) * 100 + 1,
                    "length": 100,
                }
            )

    add_variants("small", "DEL", 2)
    add_variants("medium", "DEL", 40)
    add_variants("large", "DEL", 400)
    first_small_inv_index = len(rows)
    add_variants("small", "INV", 8)
    add_variants("medium", "INV", 20)
    add_variants("large", "INV", 200)
    rows.append(dict(rows[first_small_inv_index]))  # duplicate cluster-derived row

    sv_block_path = tmp_path / "sv_block_summary.chr1.tsv"
    haploblock_path = tmp_path / "haploblocks.chr1.tsv"
    pd.DataFrame(rows).to_csv(sv_block_path, sep="\t", index=False)
    haploblocks.to_csv(haploblock_path, sep="\t", index=False)

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "genome_build": "GRCh38",
                "thresholds": {"association_probability": 0.75},
                "paths": {
                    "sv_block_summary": {"chr1": str(sv_block_path)},
                    "haploblocks": {"chr1": str(haploblock_path)},
                },
            },
            sort_keys=False,
        )
    )

    out_dir = tmp_path / "stage5_output"
    stage5_type_enrichment.main(["--config", str(config_path), "--out-dir", str(out_dir)])

    result = pd.read_csv(out_dir / "sv_type_enrichment.tsv", sep="\t")
    flagged = result[result["flagged"]]
    assert list(
        flagged[["haploblock_id", "sv_type"]].itertuples(index=False, name=None)
    ) == [("small", "INV")]

    spike = result[
        (result["haploblock_id"] == "small") & (result["sv_type"] == "INV")
    ].iloc[0]
    assert spike["observed_count"] == 8
    assert spike["expected_count"] < 2
    assert spike["q_value"] < 0.05

    totals = result.groupby("sv_type")[["observed_count", "expected_count"]].sum()
    assert np.allclose(totals["observed_count"], totals["expected_count"])

    output_config = yaml.safe_load((out_dir / "config.yaml").read_text())
    assert output_config["paths"]["sv_type_enrichment"] == str(
        (out_dir / "sv_type_enrichment.tsv").resolve()
    )


def test_stage5_optionally_collapses_same_locus_with_similar_lengths() -> None:
    haploblocks = pd.DataFrame([
        {"haploblock_id": "block", "chrom": "chr1", "start": 0, "end": 1_000},
    ])
    sv_blocks = pd.DataFrame([
        {"sv_record_id": "first", "haploblock_id": "block", "chrom": "chr1",
         "start": 100, "end": 101, "sv_type": "INS", "length": 300},
        {"sv_record_id": "similar", "haploblock_id": "block", "chrom": "chr1",
         "start": 100, "end": 101, "sv_type": "INS", "length": 305},
        {"sv_record_id": "different", "haploblock_id": "block", "chrom": "chr1",
         "start": 100, "end": 101, "sv_type": "INS", "length": 320},
    ])

    preserved = stage5_type_enrichment.enrichment_table(
        sv_blocks, haploblocks, q_threshold=0.05
    )
    collapsed = stage5_type_enrichment.enrichment_table(
        sv_blocks, haploblocks, q_threshold=0.05, collapse_length_tolerance=10
    )

    assert preserved.loc[0, "observed_count"] == 3
    assert collapsed.loc[0, "observed_count"] == 2


def test_stage5_handles_an_empty_assignment_table() -> None:
    result = stage5_type_enrichment.enrichment_table(
        pd.DataFrame(columns=["sv_id", "haploblock_id", "sv_type"]),
        pd.DataFrame([{"haploblock_id": "block1", "chrom": "chr1", "start": 0, "end": 100}]),
        q_threshold=0.05,
    )

    assert list(result.columns) == stage5_type_enrichment.ENRICHMENT_COLUMNS
    assert result.empty
