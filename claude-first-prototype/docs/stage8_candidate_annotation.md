# Stage 8: consequence-aware candidate annotation

## Purpose and inputs

Stage 8 adds genomic context to Stage 6 association candidates and Stage 7 representation candidates. It reads `sv_cluster_summary`, optional Stage 4 `sv_classification`, optional Stage 7 `sv_hash_representation`, and the GTF path registered by Stage 1. `--gtf` can replace the registered annotation.

The stage keeps association evidence, representation evidence, call quality, and gene context in separate columns. It does not calculate a single biological-importance score.

## Annotation method

The GTF is converted to 0-based, half-open gene and exon intervals. Annotation depends on SV type:

| SV type | Priority consequences |
|---|---|
| DEL | `complete_exon_loss`, `exonic_deletion`, then `intragenic_deletion`. |
| DUP | `complete_gene_duplication`, then `partial_gene_duplication`. |
| INV | `inversion_breakpoint_disruption` if either breakpoint lies in a gene; otherwise `gene_within_inversion`. |
| INS | `exonic_insertion` or `intragenic_insertion` at the insertion point. |
| Other | `gene_overlap` where relevant. |

No annotated overlap yields `intergenic_or_unannotated`. Inversions are deliberately handled differently: a gene within an inverted span is not assumed to be disrupted unless a breakpoint intersects it.

## Outputs

| File | Row grain | Important fields |
|---|---|---|
| `annotated_sv_candidates.tsv` | Candidate SV–block pair | Stage 6 summary fields, Stage 4/7 context where available, annotation fields, and call-quality label. |
| `config.yaml` | Run | Carried-forward paths plus `annotated_sv_candidates`. |

| Field | Meaning |
|---|---|
| `consequence` | Type-aware annotation label. |
| `genes` | Semicolon-delimited gene names associated with the selected annotation. |
| `overlap_basis` | Whether the result comes from a gene/exon span, an inversion breakpoint, or neither. |
| `call_quality` | `pass_precise`, `pass_imprecise`, `nonpass_precise`, or `nonpass_imprecise`, derived from the VCF fields. |
| `consequence_priority` | Sort order for consequence labels; not a severity or functional-impact score. |
| `representation_pattern`, `population_context_pattern` | Stage 7 reason the locus is a representation candidate. |

## Running and caveat

```bash
python claude-first-prototype/pipeline/stage8_candidate_annotation.py \
  --config claude-first-prototype/stage7_output/config.yaml
```

Gene and exon overlap is positional annotation, not proof of molecular or phenotypic impact. Inspect the SV representation and underlying reads before interpreting a candidate as functional.
