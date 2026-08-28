"""Stage 5: test SV-type enrichment within haploblocks.

Reads Stage 1's chromosome-specific ``sv_block_summary`` and ``haploblocks``
tables. Each unique SV-haploblock pair is counted once, independently of how
many clusters passed Stage 1's association threshold.

For each SV type, the expected count in a block is its length multiplied by
the genome-wide rate across haploblocks. Optionally, records sharing exact
coordinates, type, and sufficiently similar SV lengths are collapsed before
counting. Two-sided Poisson p-values are adjusted across the complete
haploblock-by-type grid using Benjamini-Hochberg.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("stage5_type_enrichment")

ENRICHMENT_COLUMNS = [
    "haploblock_id", "sv_type", "observed_count", "expected_count",
    "p_value", "q_value", "flagged",
]


def collapse_similar_locus_records(
    sv_blocks: pd.DataFrame,
    length_tolerance: int | None,
) -> pd.DataFrame:
    """Collapse records only when their locus, type, and length agree closely."""
    if length_tolerance is None:
        return sv_blocks
    keys = ["haploblock_id", "chrom", "start", "end", "sv_type"]
    collapsed = []
    for _, group in sv_blocks.groupby(keys, sort=False):
        group = group.sort_values("length")
        group_start = None
        for _, row in group.iterrows():
            if group_start is None or row["length"] - group_start > length_tolerance:
                collapsed.append(row)
                group_start = row["length"]
    return pd.DataFrame(collapsed, columns=sv_blocks.columns)


def resolve_path(path: str, config_dir: Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else config_dir / path


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def enrichment_table(
    sv_blocks: pd.DataFrame,
    haploblocks: pd.DataFrame,
    q_threshold: float,
    collapse_length_tolerance: int | None = None,
) -> pd.DataFrame:
    """Return the complete haploblock-by-SV-type enrichment table."""
    assigned = sv_blocks.drop_duplicates(["sv_record_id", "haploblock_id"])
    assigned = collapse_similar_locus_records(assigned, collapse_length_tolerance)
    blocks = haploblocks.drop_duplicates("haploblock_id").copy()
    blocks["length"] = blocks["end"] - blocks["start"]

    sv_types = sorted(assigned["sv_type"].unique())
    if not sv_types:
        return pd.DataFrame(columns=ENRICHMENT_COLUMNS)
    observed = (
        assigned.groupby(["haploblock_id", "sv_type"]).size()
        .unstack("sv_type", fill_value=0)
        .reindex(index=blocks["haploblock_id"], columns=sv_types, fill_value=0)
    )

    total_length = blocks["length"].sum()
    rate_by_type = observed.sum(axis=0) / total_length
    block_lengths = blocks.set_index("haploblock_id")["length"].reindex(observed.index)
    expected = np.outer(block_lengths.to_numpy(), rate_by_type.to_numpy())

    observed_flat = observed.to_numpy().ravel().astype(float)
    expected_flat = expected.ravel()
    lower_tail = stats.poisson.cdf(observed_flat, expected_flat)
    upper_tail = stats.poisson.sf(observed_flat - 1, expected_flat)
    p_values = np.minimum(1.0, 2 * np.minimum(lower_tail, upper_tail))
    q_values = benjamini_hochberg(p_values)

    rows = []
    for block_index, haploblock_id in enumerate(observed.index):
        for type_index, sv_type in enumerate(sv_types):
            flat_index = block_index * len(sv_types) + type_index
            rows.append(
                {
                    "haploblock_id": haploblock_id,
                    "sv_type": sv_type,
                    "observed_count": int(observed_flat[flat_index]),
                    "expected_count": float(expected_flat[flat_index]),
                    "p_value": float(p_values[flat_index]),
                    "q_value": float(q_values[flat_index]),
                    "flagged": bool(q_values[flat_index] < q_threshold),
                }
            )

    return pd.DataFrame(rows, columns=ENRICHMENT_COLUMNS).sort_values(
        ["q_value", "haploblock_id", "sv_type"]
    ).reset_index(drop=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("stage1_output/config.yaml"))
    parser.add_argument("--out-dir", type=Path, default=Path("stage5_output"))
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument(
        "--collapse-length-tolerance", type=int, default=None,
        help=(
            "Collapse records with the same chrom/start/end/type when their "
            "length differs by at most this many bp; default preserves every record."
        ),
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text())
    config_dir = args.config.parent

    sv_blocks = pd.concat(
        [
            pd.read_csv(resolve_path(path, config_dir), sep="\t")
            for path in config["paths"]["sv_block_summary"].values()
        ],
        ignore_index=True,
    )
    haploblocks = pd.concat(
        [
            pd.read_csv(resolve_path(path, config_dir), sep="\t")
            for path in config["paths"]["haploblocks"].values()
        ],
        ignore_index=True,
    )
    result = enrichment_table(
        sv_blocks, haploblocks, args.q_threshold, args.collapse_length_tolerance
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.out_dir / "sv_type_enrichment.tsv"
    result.to_csv(output_path, sep="\t", index=False)
    log.info(
        "Tested %d haploblock-by-type cells; %d flagged at q < %.3g",
        len(result),
        int(result["flagged"].sum()),
        args.q_threshold,
    )

    stage5_config = dict(config)
    stage5_config["thresholds"] = dict(config["thresholds"])
    stage5_config["thresholds"]["sv_type_enrichment_q"] = args.q_threshold
    stage5_config["settings"] = dict(config.get("settings", {}))
    stage5_config["settings"]["collapse_length_tolerance_bp"] = (
        args.collapse_length_tolerance
    )
    stage5_config["paths"] = dict(config["paths"])
    stage5_config["paths"]["sv_type_enrichment"] = str(output_path.resolve())
    (args.out_dir / "config.yaml").write_text(yaml.safe_dump(stage5_config, sort_keys=False))


if __name__ == "__main__":
    main(sys.argv[1:])
