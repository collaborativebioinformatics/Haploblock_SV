# Changelog

## Unreleased

### Overview

- Reconciled the early pipeline around stable Stage 1 contracts for SV genotypes, sample metadata, haploblocks, cluster memberships, unique SV–block summaries, and passing SV–cluster associations.
- Clarified Stage 2 as an optional descriptive boundary analysis, separate from cluster association.
- Ported Stage 4 to classify common and population-specific SVs directly from Stage 1 genotypes and population metadata, without using haploblock cluster labels.
- Ported Stage 5 to test per-haploblock SV-type enrichment from unique SV–block pairs using block-length-adjusted Poisson expectations and Benjamini–Hochberg correction.
- Added focused contract tests and carried-forward configuration files so later stages can consume the reconciled outputs consistently.

### Details

- Restored indexed Stage 2 interval lookup, validated non-overlapping haploblocks, and report proximity for SVs just outside a block.
- Kept Stage 2 outputs schema-stable when a requested chromosome has no SV records.
- Reduced Stage 4 memory use by aggregating genotype counts one sample at a time; multiallelic genotype alleles are now parsed correctly.
- Split Stage 4's per-population and per-SV output construction so neither requires deduplication of the other.
- Made Stage 5 emit a schema-correct empty enrichment table when there are no assigned SVs.
- Replaced Stage 1's argv reconstruction with a typed pipeline call, extracted shared SV schema helpers from the executable module, removed unused per-chromosome QC return values, and relaxed the Stage 1 contract test to allow future published paths.
