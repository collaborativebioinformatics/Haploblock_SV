"""Generate a tiny Stage-2-style example for pipeline/stage5_type_enrichment.py.

Deterministic (fixed seed). Six contiguous haploblocks of very different
lengths; DEL and INS counts are drawn Poisson-proportional to block length
(so no block should be flagged for those), and a small block gets an
artificial INV spike (which should be flagged after FDR).

Run:  python example_data/stage5_example/make_example.py
Writes sv_calls.tsv and haploblocks.tsv next to this file.
"""

from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(7)

# haploblock_id, start, end  (0-based half-open, contiguous)
BLOCKS = [
    ("hb1", 0, 12_000),
    ("hb2", 12_000, 32_000),
    ("hb3", 32_000, 62_000),
    ("hb4", 62_000, 110_000),
    ("hb5", 110_000, 180_000),
    ("hb6_small", 180_000, 180_800),   # tiny block -> carries the INV spike
]

hb = pd.DataFrame(
    [
        {"haploblock_id": b, "chrom": "chr1", "start": s, "end": e,
         "n_snps": np.nan, "n_clusters": np.nan, "hash_length": np.nan,
         "cluster_diff_score": np.nan}
        for b, s, e in BLOCKS
    ]
)

# background rates (per bp) -- DEL/INS counts come out proportional to length
RATES = {"DEL": 1 / 1_200, "INS": 1 / 2_500, "INV": 1 / 60_000}

rows = []
sv_i = 0


def add(start, end, sv_type, hbid, pos="within_block"):
    global sv_i
    rows.append({
        "sv_id": f"sv{sv_i:04d}", "chrom": "chr1", "start": int(start), "end": int(end),
        "sv_type": sv_type, "imprecise": False,
        "length": 300 if sv_type == "INS" else int(end - start),
        "S1": "0/1", "S2": "0/0",                    # dummy genotype cols (Stage 5 ignores them)
        "position_class": pos, "haploblock_id": hbid,
    })
    sv_i += 1


for b, s, e in BLOCKS:
    length = e - s
    for sv_type, rate in RATES.items():
        for _ in range(rng.poisson(rate * length)):
            p = int(rng.integers(s, max(s + 1, e - 1)))
            end = p + 1 if sv_type == "INS" else p + int(rng.integers(50, 400))
            add(p, min(end, e), sv_type, b)

# artificial INV spike in the small block
for _ in range(6):
    p = int(rng.integers(180_000, 180_780))
    add(p, p + 100, "INV", "hb6_small")

# one boundary-crossing SV (comma-joined id -> Stage 5 explodes it into both blocks)
add(31_950, 32_050, "DEL", "hb2,hb3", pos="boundary_crossing")
# one outside_block SV (no haploblock -> Stage 5 drops it)
add(500_000, 500_200, "DEL", "", pos="outside_block")

sv = pd.DataFrame(rows)
sv.to_csv(HERE / "sv_calls.tsv", sep="\t", index=False)
hb.to_csv(HERE / "haploblocks.tsv", sep="\t", index=False)
print(f"wrote {len(sv)} SVs across {len(hb)} haploblocks to {HERE}")
print(sv.groupby(["haploblock_id", "sv_type"]).size().unstack(fill_value=0))
