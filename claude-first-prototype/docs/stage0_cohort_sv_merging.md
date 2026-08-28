# Stage 0: cohort SV merging

## Purpose

Stage 0 will turn single-sample long-read SV calls into a cohort VCF in which records are intended to represent comparable biological alleles. This is the prerequisite for interpreting an SV on multiple haploblock backgrounds: without record reconciliation, that pattern could reflect different calls at a repeat rather than one allele.

## Status and inputs

This stage is not implemented in the prototype. Stages 1–9 currently use the pre-merged cohort VCF at `input/1kgp_ont_cohort.postfilter.full.vcf.gz`.

The planned input is a list of single-sample long-read SV VCFs. Each VCF must use the same genome build and have a stable sample identifier. Reconciliation must be done separately by SV type where the underlying comparison method requires it.

## Planned method

1. Merge the sample VCFs into a multi-sample callset.
2. Use `truvari collapse` to consolidate records judged to describe the same event.
3. Apply cohort-level record filtering and preserve the source and filter information needed to evaluate calls later.

The current input was produced with Sniffles calling, `bcftools merge`, `truvari collapse`, and kanpig regenotyping. Kanpig requires resolved variant sequences, so that particular callset contains DEL and INS records but not INV or BND records.

## Expected output contract

| VCF element | Requirement | Downstream use |
|---|---|---|
| Record identifier | Unique and stable within the VCF | Stage 1 creates its own unique `sv_record_id`, but preserves this as `sv_id`. |
| Chromosome and coordinates | GRCh38, with an SV type and usable span | Defines haploblock overlap and gene context. |
| `SVTYPE` and `SVLEN`/`END` | Present where applicable | Used for type-specific analyses and annotation. |
| `FILTER` and `IMPRECISE` | Preserved rather than silently removed | Retained as call-quality context. |
| Per-sample `GT` | Present and interpretable | Used to calculate dosage, frequency, and cluster association. |

## Major limitations

Calls in low-complexity or repetitive sequence may remain difficult to reconcile. An apparent multi-cluster SV should therefore be checked for local representation and read-level support before it is interpreted as recurrent mutation, recombination, or an old allele.
