# Stage 5: per-haploblock SV-type enrichment

## Purpose

Stage 5 is an optional screen for haploblocks containing more or fewer SVs of a particular type than expected from their length. A flagged result can motivate a hypothesis about local sequence architecture, but does not establish a mechanism or test haploblock representation.

## Inputs and calculation

The stage reads all chromosome-specific `sv_block_summary` and `haploblocks` tables from Stage 1. It counts every unique `sv_record_id, haploblock_id` pair once and creates the full grid of haploblocks and observed SV types, including zero-count cells.

For each SV type, the expected count in a block is:

`overall count of that SV type × block length / total haploblock length`

A two-sided Poisson probability measures how unusual the observed count is under this length-only model. Benjamini–Hochberg adjustment produces a q-value across every tested block/type cell. In plain language, the q-value controls the expected proportion of false-positive flags among the cells called significant.

## Outputs

| File | Row grain | Fields |
|---|---|---|
| `sv_type_enrichment.tsv` | Haploblock–SV-type combination | `haploblock_id`, `sv_type`, `observed_count`, `expected_count`, `p_value`, `q_value`, `flagged`. |
| `config.yaml` | Run | Carried-forward paths, q-value threshold, and any record-collapse tolerance. |

`flagged` is true when `q_value < --q-threshold` (default 0.05). A p-value is the result for one block/type cell; a q-value reflects that many cells were screened.

## Running and caveats

```bash
python claude-first-prototype/pipeline/stage5_type_enrichment.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

By default every VCF record remains distinct. `--collapse-length-tolerance BP` collapses only records with the same chromosome, start, end, and SV type whose lengths differ by no more than `BP`. The model adjusts only for block length. Differences in callability, local repeat content, SV count variability, and record representation can all affect the result.
