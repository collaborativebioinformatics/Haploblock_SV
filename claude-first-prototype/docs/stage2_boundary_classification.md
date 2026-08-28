# Stage 2: boundary classification

## Purpose and inputs

Stage 2 is an optional spatial description of SVs relative to haploblock boundaries. It reads `haploblocks` and, when available, `sv_genotypes` from the Stage 1 configuration. It does not use cluster association calls and therefore cannot establish whether a haploblock predicts an SV.

## Method

For every SV interval, the stage counts exactly overlapping haploblocks. It also tests whether either endpoint lies within the configured distance of any block boundary. The default distance is inherited from Stage 1 (`5000` bp) unless `--boundary-distance-bp` overrides it.

| `position_class` | Rule |
|---|---|
| `outside_block` | Overlaps no haploblock. |
| `within_block` | Overlaps exactly one haploblock. |
| `boundary_crossing` | Overlaps two or more haploblocks. |

`near_boundary` is separate from this classification. A nearby SV is not relabeled as boundary crossing.

## Outputs

| File | Row grain | Important fields |
|---|---|---|
| `boundary_svs.<chrom>.tsv` | SV record | Fixed Stage 1 SV metadata plus `position_class`, `haploblock_id`, and `near_boundary`. Multiple overlapping IDs are comma-separated. |
| `boundary_qc.json` | Run and chromosome | Boundary-distance setting, number of SVs, counts by `position_class`, near-boundary count, and output path. |

## Running and caveat

```bash
python claude-first-prototype/pipeline/stage2_intersect.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

An SV that crosses a boundary may identify a locus that deserves follow-up, but its position alone does not show that the SV is poorly represented by haploblock clusters.
