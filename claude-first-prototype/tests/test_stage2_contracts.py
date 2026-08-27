"""Contract tests for exact crossing and boundary proximity."""

import sys
from pathlib import Path

import pandas as pd

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

from stage2_intersect import classify_sv_positions


def test_crossing_and_proximity_are_separate() -> None:
    haploblocks = pd.DataFrame(
        [
            {"haploblock_id": "block1", "chrom": "chr1", "start": 0, "end": 100},
            {"haploblock_id": "block2", "chrom": "chr1", "start": 100, "end": 200},
        ]
    )
    variants = pd.DataFrame(
        [
            {"sv_id": "inside", "chrom": "chr1", "start": 20, "end": 30},
            {"sv_id": "near", "chrom": "chr1", "start": 95, "end": 99},
            {"sv_id": "crossing", "chrom": "chr1", "start": 95, "end": 105},
            {"sv_id": "outside", "chrom": "chr1", "start": 250, "end": 260},
        ]
    )

    result = classify_sv_positions(variants, haploblocks, boundary_bp=10).set_index("sv_id")

    assert result.loc["inside", "position_class"] == "within_block"
    assert not result.loc["inside", "near_boundary"]
    assert result.loc["near", "position_class"] == "within_block"
    assert result.loc["near", "near_boundary"]
    assert result.loc["crossing", "position_class"] == "boundary_crossing"
    assert result.loc["crossing", "haploblock_id"] == "block1,block2"
    assert result.loc["crossing", "near_boundary"]
    assert result.loc["outside", "position_class"] == "outside_block"
    assert not result.loc["outside", "near_boundary"]
