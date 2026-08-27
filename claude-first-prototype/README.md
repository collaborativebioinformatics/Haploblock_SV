# Haploblock_SV — Structural Variants Within Haploblocks

Prototype pipeline plan for the Structural Variants Hackathon at Baylor College of Medicine, August 25–28, 2026. This document preserves the original brainstorming and planning pass (see `PROMPTS.md`) while identifying which early stages are currently implemented. Later analyses remain proposals while the biological plan evolves.

## Introduction

Haploblocks — genomic regions of conserved haplotype structure identified by the [haploblocks.org](https://haploblocks.org) / [data.haploblocks.org](https://data.haploblocks.org) projects via genomic hashing — capture a layer of population structure that is complementary to single-SNP analyses, but how structural variants (SVs) are organized relative to these blocks is unexplored. Our core research question is whether SV type, location, and population specificity are non-randomly distributed across haploblocks, and whether SV-bearing haplotypes correspond to the clusters already derived from small-variant haplotype hashes. The current primary input is a pre-merged 1000 Genomes ONT VCF; proposed enrichment, spatial-statistics, PCA/UMAP, cluster-agreement, and functional-overlay analyses below are retained as possible follow-up work rather than current commitments. What's novel is treating haploblock **clusters**, rather than only immutable block regions, as the unit of SV interpretation.

## Pipeline overview

| # | Stage | Purpose | Addresses |
|---|---|---|---|
| 0 | Cohort SV merging | Download/accept single-sample long-read VCFs, merge samples, and reconcile equivalent representations with `truvari collapse` | Placeholder; Linh is working on it |
| 1 | Cluster-aware preprocessing | Infer SV-to-haploblock-cluster associations from a merged cohort VCF and write probability/support/QC outputs | Implemented |
| 2 | Boundary classification | Descriptively classify each SV as safely within or near/crossing a block boundary | Optional/implemented |
| 3 | Boundary enrichment test | Earlier proposal for permutation/spacing tests; not part of the current pipeline | Future decision |
| 4 | Common vs. population-specific SV classification | Calculate AF using populations supplied by `sample_metadata.tsv`, independently of haploblocks.org clusters | Proposed |
| 5 | Per-haploblock SV-type enrichment | Poisson/negative-binomial regression, block length + SNP density as offset, per-block-per-type deviation, BH-FDR corrected, minimum-count threshold | H2, H4 |
| 6 | Population-cluster correlation | Compare population-specific SV patterns with predefined SNV-based haplotype clusters | Proposed |
| 7 | SV-based population structure reconstruction | Per-sample/per-haploblock SV matrix → PCA/UMAP and cluster-agreement summaries | Proposed visualization |
| 8 | Duplication/inversion gene overlay | Gene overlap for supported recurrent/population-specific SV types; INV depends on a future Stage 0 that retains inversions | Proposed |
| 9 | Integration & report | Aggregate analysis results into a per-haploblock summary and plots | Proposed |

If pursued, Stages 4–8 contain separable questions and can largely be divided across the team after Stage 1 outputs have been evaluated.

## Original hackathon development plan (historical)

The schedule below is preserved from the initial planning pass. It is not the current execution plan; in particular, current Stage 0/1 are defined in the implementation sections below and complex statistical/visualization stages remain future decisions.

**Day 1 (Aug 25).** Original goal: stand up ingestion, QC, and intersection. Current replacement: Linh is developing Stage 0 merging, Stage 1 performs cluster-aware preprocessing from the existing pre-merged VCF, and Stage 2 is optional boundary classification.

**Day 2 (Aug 26 — today).** Stages 3–5: boundary enrichment test, population-specific classification, and the size-adjusted per-block SV-type enrichment scan. This is the statistical core of the project — get p-values/FDR-corrected results on the real dataset today, even if plots/report are rough.

**Day 3 (Aug 27).** Stages 6–8 in parallel across the team (population-cluster correlation, SV-based clustering + ARI, gene overlay for Maria/Alistair's duplication/inversion question). Start Stage 9 integration in the afternoon so there's a working summary table before the final day.

**Day 4 / buffer (Aug 28).** Original goal: finish Stage 9, run an end-to-end integration check, and prepare the demo/writeup. The previously suggested web view is outside the current scope.

## Proposed testing and validation ideas

- **Per-stage sanity checks on a small slice:** once the downstream analyses are selected, run those stages on one chromosome (e.g., chr22) or a synthetic haploblock+SV set before running genome-wide; confirm intersection counts (Stage 2), classification counts (Stage 4), and enrichment p-value distributions (Stage 5) look sane (e.g., p-values roughly uniform under a shuffled-label negative control).
- **Negative controls:** re-run Stage 3's permutation test and Stage 5's regression on label-shuffled data — both should show no significant enrichment; a positive result on shuffled data indicates a pipeline bug, not a biological signal.
- **End-to-end integration run:** once Stage 0 is implemented, execute the selected stages on one chromosome and confirm the expected outputs are populated.
- **Validation against known structure:** if Stage 7 is pursued, compare SV-based structure with population labels supplied in `sample_metadata.tsv` before interpreting agreement with haploblocks.org clusters.
- **Reproducibility check:** re-run the full pipeline twice with the same config (fixed seeds) and confirm identical output; re-run once with a different seed and confirm permutation/UMAP results are stable within expected tolerance.

## Non-pip dependencies

- `bcftools` and `truvari` — planned Stage 0 merge/collapse workflow being developed by Linh.
- R is not required; the plan above is Python-only. If a team member prefers R for Stage 5/6 stats (e.g. `MASS::glm.nb`), that's a drop-in alternative to `statsmodels`, not a pipeline dependency.

## Current Stage 0 placeholder: cohort SV merging

Stage 0 is not yet implemented in this prototype. **Linh is working on it.** Its eventual role is to download or accept a list of single-sample long-read SV VCFs, combine the samples, and reconcile calls that describe the same biological SV differently. This representation-merging step is necessary because basecalling errors, mapping ambiguity around repeats, and variant-caller differences can shift breakpoints or otherwise produce different records for the same event.

The planned Stage 0 will use sample merging followed by `truvari collapse`. It should also permit a different cohort to be supplied as a list of single-sample VCFs. It will not run kanpig. Consequently, it may retain inversions, while BND records will probably remain outside the merged callset.

For the hackathon, Stage 1 starts from the existing pre-merged VCF at `input/1kgp_ont_cohort.postfilter.full.vcf.gz`. That file was produced previously using:

1. Sniffles SV calling
2. `bcftools merge`
3. `truvari collapse`
4. kanpig regenotyping after collapsing equivalent representations

Kanpig needs resolved variant sequences, so this current input contains DEL and INS records, but no INV or BND records. This limitation describes the current input, not the intended general scope of future analyses.

## Current Stage 1: cluster-aware preprocessing

`pipeline/stage1_cluster_aware.py` is the implemented entry point. It reads the merged VCF directly, downloads or reuses the per-block cluster membership files from data.haploblocks.org, and infers which haploblock cluster or clusters carry each SV.

```bash
python claude-first-prototype/pipeline/stage1_cluster_aware.py
```

The default output directory is `claude-first-prototype/stage1_output/`. Chromosomes are processed separately so work can be parallelized and rerun selectively. Stage 1 now keeps the normalized inputs needed by later analyses instead of treating them as temporary files:

- `sv_genotypes.<chrom>.tsv`: one row per SV, with normalized metadata followed by one raw GT column per sample. This is the Stage 4 and Stage 7 genotype contract.
- `samples.tsv`: canonical sample IDs and their original VCF IDs.
- `haploblocks.<chrom>.tsv`: one row per haploblock.
- `cluster_memberships.<chrom>.tsv`: one row per represented sample haplotype and haploblock cluster. This is the independent cluster input for Stages 6 and 7.
- `sv_block_summary.<chrom>.tsv`: one row per overlapping SV and haploblock, including association counts but never duplicating a pair by passing cluster. This is the counting input for Stage 5.
- `sv_to_clusters.<chrom>.tsv`: one row per SV, haploblock, and cluster that passes the association threshold.

All paths are registered in `stage1_output/config.yaml`. Population labels remain independent of cluster inference: pass the Stage 0 table with `--sample-metadata` to register it for Stages 4, 6, and 7. Useful method diagnostics are kept under `debug_and_qc/`; downloaded cluster files remain temporary and are removed after a successful run.

The table keys and row meanings are part of the contract:

| Path key | One row per | Required identity columns |
|---|---|---|
| `sv_genotypes` | input VCF record | `sv_id`, `chrom`, `start`, `end`, `sv_type`; sample GT columns follow the fixed metadata columns |
| `haploblocks` | haploblock | `haploblock_id`, `chrom`, `start`, `end` |
| `cluster_memberships` | complete sample haplotype assignment | `haploblock_id`, `sample_id`, `haplotype`, `cluster_id` |
| `sv_block_summary` | overlapping SV–haploblock pair | `sv_id`, `haploblock_id`; this key is unique even when several clusters pass |
| `sv_to_clusters` | passing SV–haploblock–cluster association | `sv_id`, `haploblock_id`, `cluster_id` |

The optional `sample_metadata` table must use the canonical `sample_id` values in `samples.tsv` and provide a `population` column. A `superpopulation` column may also be supplied, but downstream analyses must state explicitly which grouping they use.

### Downstream ownership

| Stage | Stage 1 inputs |
|---|---|
| 4 | `sv_genotypes`, `sample_metadata`, and `sv_block_summary` |
| 5 | `sv_block_summary` and `haploblocks`; count unique `sv_id, haploblock_id` pairs |
| 6 | Stage 4 classifications plus `cluster_memberships` and `sample_metadata` |
| 7 | `sv_genotypes`, `cluster_memberships`, and `sample_metadata` |
| 8 | Stage 4 classifications plus SV coordinates and types from `sv_genotypes` |

### Interpreting `sv_to_clusters.tsv`

Each row is an **SV–haploblock–cluster association**, not merely an overlap between an SV and an immutable block region. An SV can therefore have no passing cluster, one passing cluster, multiple passing clusters in the same block, or associations in more than one overlapped block. Absence from this table means that no evaluated cluster passed the configured probability and support thresholds; it does not prove that the SV is absent.

The main fields are:

| Field | Interpretation |
|---|---|
| `sv_id` | Input VCF ID when it is short enough. IDs longer than the configured limit are replaced with `SV_<chrom>_<start>_<end>_<type>_<hash>`, where the hash preserves uniqueness. |
| `chrom`, `start`, `end` | SV coordinates in 0-based, half-open convention. |
| `sv_type`, `length`, `filter`, `imprecise` | Metadata retained from or derived from the input VCF. |
| `block_start`, `block_end` | Haploblock region whose cluster membership was evaluated. These coordinates identify the context; they are not themselves evidence of association. |
| `cluster_id` | Representative haplotype string that identifies the cluster. |
| `cluster_haplotypes_total` | Total haplotypes listed in the cluster membership file. |
| `cluster_haplotypes_in_vcf` | Cluster haplotypes whose canonicalized sample is represented in the VCF. |
| `callable_haplotypes` | Number of represented cluster haplotypes with non-missing VCF evidence. Genotypes containing `.` are excluded rather than treated as reference. |
| `required_callable_haplotypes` | Adaptive support requirement: the smaller of three or the number of VCF-represented haplotypes in the cluster. |
| `call_rate` | Fraction of VCF-represented cluster haplotypes that contributed callable evidence. |
| `expected_alt_haplotypes` | Expected number of alternate-bearing haplotypes, including probabilistic assignments for heterozygotes. |
| `sv_probability` | EM estimate of the probability that a haplotype in this cluster carries the SV. Values at or above the association threshold (default 0.75) are emitted here. |
| `ci95_low`, `ci95_high` | Approximate 95% interval around the estimated cluster probability. |
| `evidence_tier` | `low` when fewer than three callable haplotypes support the estimate; otherwise `standard`. Low-evidence rows can be biologically interesting but should be interpreted cautiously. |
| `model_converged`, `model_iterations` | Whether the EM calculation met its tolerance and how many iterations it used. |
| `overlaps_multiple_blocks` | Indicates that the SV overlaps more than one haploblock and was evaluated independently in each one. |

For heterozygous samples, the VCF phase is not assumed to match the cluster-generation `hap0`/`hap1` labels. Assignment instead uses the other callable members of both clusters. A heterozygote is assigned to one cluster only when its posterior probability reaches the threshold; otherwise it remains ambiguous. Cluster `hap0` and `hap1` labels do not imply maternal and paternal origin.

## Optional Stage 2: boundary classification

`pipeline/stage2_intersect.py` remains available for the separate descriptive question of whether an SV lies safely within a block or near/crosses a block boundary. It reads Stage 1's generated config:

```bash
python claude-first-prototype/pipeline/stage2_intersect.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

This analysis is not required to establish cluster association. Because published blocks are contiguous regions, boundary results should be treated as a distinct spatial analysis rather than evidence that an SV belongs to a cluster.

## Current testing status

The earlier dbVar ingestion/QC tests described in the historical prototype no longer match the current pipeline. Tests for the cluster-aware workflow will be handled in a separate testing pass.
