# Haploblock_SV — Structural Variants Within Haploblocks

Prototype pipeline plan for the Structural Variants Hackathon at Baylor College of Medicine, August 25–28, 2026. This document preserves the original brainstorming and planning pass (see `PROMPTS.md`) while identifying which early stages are currently implemented. Later analyses remain proposals while the biological plan evolves.

## Introduction

Haploblocks — genomic regions of conserved haplotype structure identified by the [haploblocks.org](https://haploblocks.org) / [data.haploblocks.org](https://data.haploblocks.org) projects via genomic hashing — capture a layer of population structure that is complementary to single-SNP analyses, but how structural variants (SVs) are organized relative to these blocks is unexplored. Our core research question is narrower than asking whether SVs recapitulate ancestry: **what information about structural variation is captured, or missed, by these particular recombination-defined genomic-hash clusters?** The current primary input is a pre-merged 1000 Genomes ONT VCF; spatial-statistics, PCA/UMAP, cluster-agreement, and functional-overlay analyses below are retained as possible follow-up work rather than current commitments. What's novel is treating haploblock **clusters**, rather than only immutable block regions, as the unit of SV interpretation and directly measuring when those clusters do or do not resolve SV carriage.

## Pipeline overview

| # | Stage | Purpose | Addresses |
|---|---|---|---|
| 0 | Cohort SV merging | Download/accept single-sample long-read VCFs, merge samples, and reconcile equivalent representations with `truvari collapse` | Placeholder; Linh is working on it |
| 1 | Cluster-aware preprocessing | Infer SV-to-haploblock-cluster associations from a merged cohort VCF and write probability/support/QC outputs | Implemented |
| 2 | Boundary classification | Descriptively classify each SV as safely within or near/crossing a block boundary | Optional/implemented |
| 3 | Boundary enrichment test | Earlier proposal for permutation/spacing tests; not part of the current pipeline | Future decision |
| 4 | Common vs. population-specific SV classification | Calculate AF using populations supplied by `sample_metadata.tsv`, independently of haploblocks.org clusters | Implemented |
| 5 | Per-haploblock SV-type enrichment | Length-adjusted Poisson tests across the complete block-by-type grid, with BH-FDR correction | Implemented |
| 6 | Population-conditioned SV–cluster association | Test whether local SNV-derived clusters predict SV carriage beyond population membership and whether associations transfer across populations | Implemented |
| 7 | Hash representation audit and structure QC | Measure carrier concentration across clusters and carrier purity within clusters; retain information gain and SV PCA as secondary diagnostics | Implemented |
| 8 | Consequence-aware representation annotation | Add type-aware gene/exon and breakpoint context to SVs captured or missed by the hashes | Implemented |
| 9 | Integration & report | Aggregate analysis results into a per-haploblock summary and plots | Proposed |

If pursued, Stages 4–8 contain separable questions and can largely be divided across the team after Stage 1 outputs have been evaluated.

### Rationale for the Stage 6–8 refinement

The initial downstream plan emphasized agreement metrics and population-structure visualizations.
Those remain useful validation and QC, but population structure from SV genotypes and population
differentiation of SVs are already expected from established population genetics. The refined plan
therefore focuses on the contribution specific to this project: determining when the local
SNV-derived genomic hashes provide portable tags for long-read SVs, when an association depends on
population background, and when SVs subdivide an existing hash and add previously unrepresented
local information. Candidate annotation then distinguishes plausible breakpoint, exon, dosage, and
span consequences rather than treating every gene overlap equally. These changes retain the useful
components of the original plan while producing more directly interpretable locus-level results.

## Biological decisions the pipeline should support

The pipeline is an evaluation of a representation, not a metric-generation exercise. Its central
outputs should support the following decisions:

| Biological question | Current coverage | What is still needed | Why we care |
|---|---|---|---|
| Are SVs reliably tagged by the same genomic-hash cluster across populations? | Stage 6 estimates population-conditioned association and directional consistency. | Use batch-aware maximum-statistic permutations and held-out or leave-one-population-out validation; distinguish positive carrier tags from negative exclusion tags. | A portable tag would let the existing hash stand in for an expensive SV genotype in cohorts where long-read data are unavailable. |
| Do SVs subdivide an apparently homogeneous SNV-derived cluster? | Stage 7 reports count-supported mixed carrier/non-carrier diplotype candidates. | Confirm the split across populations or batches and calibrate it against a within-diplotype null before describing it as stable. | A supported subdivision candidate identifies information potentially missing from the hash and tells us where direct SV measurement adds resolution. |
| Does the same SV allele occur on multiple divergent haplotype backgrounds? | Stage 1 can emit multiple cluster associations, but later stages do not classify this pattern. | First resolve multiallelic/VNTR loci, then measure the number and divergence of carrier clusters and inspect sequence/read evidence. | A validated multi-background allele can indicate recurrent mutation, gene conversion, recombination, or an old allele; an unvalidated one often exposes call collapsing or representation error. |
| Is a population-enriched SV carried on an otherwise cosmopolitan haplotype background? | Stage 7 joins Stage 4 classification to the population breadth of the top supported carrier cluster and emits a candidate flag. | Treat the flag as descriptive until adequate population and cluster support are confirmed. | This pattern is consistent with a relatively recent event that the older SNV-defined hash cannot resolve and is more informative than ancestry concordance alone. |

Two complementary measurements anchor interpretation. The **SV-centric** view asks whether carriers
of one resolved SV are concentrated in one or a few hash clusters or distributed across many clusters.
The **cluster-centric** view asks whether members of one hash cluster are homogeneous for SV carriage
or split into carriers and non-carriers. Haploblocks are sufficient SV tags where both concentration
and within-cluster purity are high. Direct SV measurement adds information where the same SV spans
many clusters or where SV presence repeatedly subdivides a cluster. Population labels are used to
control stratification and describe whether these patterns are shared or population-restricted; they
are not an alternative SV-prediction model.

This framing makes Stages 1, 4, 6, and 7 the core representation audit. Stage 8 interprets the most
credible captured and missed variants. Boundary and per-block count enrichment (Stages 2, 3, and 5)
remain optional side questions unless they produce a specific mechanistic hypothesis.

## Original hackathon development plan (historical)

The schedule below is preserved from the initial planning pass. It is not the current execution plan; in particular, current Stage 0/1 are defined in the implementation sections below and complex statistical/visualization stages remain future decisions.

**Day 1 (Aug 25).** Original goal: stand up ingestion, QC, and intersection. Current replacement: Linh is developing Stage 0 merging, Stage 1 performs cluster-aware preprocessing from the existing pre-merged VCF, and Stage 2 is optional boundary classification.

**Day 2 (Aug 26).** Stages 4–5 now provide population-specific classification and a length-adjusted per-block SV-type enrichment scan. Stage 3 boundary enrichment remains a separate future decision.

**Day 3 (Aug 27).** Stages 6–8 in parallel across the team (population-cluster correlation, SV-based clustering + ARI, gene overlay for Maria/Alistair's duplication/inversion question). Start Stage 9 integration in the afternoon so there's a working summary table before the final day.

**Day 4 / buffer (Aug 28).** Original goal: finish Stage 9, run an end-to-end integration check, and prepare the demo/writeup. The previously suggested web view is outside the current scope.

## Proposed testing and validation ideas

- **Per-stage sanity checks on a small slice:** once the downstream analyses are selected, run those stages on one chromosome (e.g., chr22) or a synthetic haploblock+SV set before running genome-wide; confirm intersection counts (Stage 2), classification counts (Stage 4), and enrichment p-value distributions (Stage 5) look sane (e.g., p-values roughly uniform under a shuffled-label negative control).
- **Negative controls:** if Stage 3 is pursued, run its permutation test on shuffled data; also rerun Stage 5 after shuffling SV types among blocks. Neither should show systematic enrichment.
- **End-to-end integration run:** once Stage 0 is implemented, execute the selected stages on one chromosome and confirm the expected outputs are populated.
- **Validation against known structure:** if Stage 7 is pursued, compare SV-based structure with population labels supplied in `sample_metadata.tsv` before interpreting agreement with haploblocks.org clusters.
- **Reproducibility check:** re-run the full pipeline twice with the same config (fixed seeds) and confirm identical output; re-run once with a different seed and confirm permutation/UMAP results are stable within expected tolerance.

## Non-pip dependencies

- `bcftools` and `truvari` — planned Stage 0 merge/collapse workflow being developed by Linh.
- R is not required; the implemented pipeline is Python-only. Stage 5 uses SciPy for its current Poisson tests.

## Current Stage 0 placeholder: cohort SV merging

Stage 0 is not yet implemented in this prototype. **Linh is working on it.** Its eventual role is to download or accept a list of single-sample long-read SV VCFs, combine the samples, and reconcile calls that describe the same biological SV differently. This representation-merging step is necessary because basecalling errors, mapping ambiguity around repeats, and variant-caller differences can shift breakpoints or otherwise produce different records for the same event.

The planned Stage 0 will use sample merging followed by `truvari collapse`. It should also permit a different cohort to be supplied as a list of single-sample VCFs. It will not run kanpig. Consequently, it may retain inversions, while BND records will probably remain outside the merged callset.

Stage 0 owns representation-level QC: reconciling equivalent calls, deduplicating them, and applying any cohort-wide FILTER or size policy. Stage 1 preserves every record it receives and reports genotype/association QC without silently removing IMPRECISE calls or particular SV types. Analysis-specific size or confidence subsets belong in the downstream analysis that requests them.

**Why we care:** an apparent SV on several haplotype backgrounds is biologically interesting only after Stage 0 has shown that it is one comparable allele rather than several repeat alleles or collapsed call representations.

For the hackathon, Stage 1 starts from the existing pre-merged VCF at `input/1kgp_ont_cohort.postfilter.full.vcf.gz`. That file was produced previously using:

1. Sniffles SV calling
2. `bcftools merge`
3. `truvari collapse`
4. kanpig regenotyping after collapsing equivalent representations

Kanpig needs resolved variant sequences, so this current input contains DEL and INS records, but no INV or BND records. This limitation describes the current input, not the intended general scope of future analyses.

## Current Stage 1: cluster-aware preprocessing

`pipeline/stage1_cluster_aware.py` is the implemented entry point. It reads the merged VCF directly, downloads or reuses the per-block cluster membership files from data.haploblocks.org, downloads the matching 1000 Genomes metadata and Ensembl GRCh38 gene annotation unless local paths are supplied, and infers which haploblock cluster or clusters carry each SV.

```bash
python claude-first-prototype/pipeline/stage1_cluster_aware.py
```

The default output directory is `claude-first-prototype/stage1_output/`. Chromosomes are processed separately so work can be parallelized and rerun selectively. Stage 1 now keeps the normalized inputs needed by later analyses instead of treating them as temporary files:

- `sv_genotypes.<chrom>.tsv`: one row per input VCF record, with normalized metadata followed by one raw GT column per sample. `sv_record_id` is the unique record key; `sv_id` is retained as a potentially repeated source label. This is the Stage 4 and Stage 7 genotype contract.
- `samples.tsv`: canonical sample IDs and their original VCF IDs.
- `sample_metadata.tsv`: a VCF-ordered copy of the automatically downloaded metadata (or `--sample-metadata`) containing represented samples and using the canonical IDs from `samples.tsv`.
- `Homo_sapiens.GRCh38.115.gtf.gz`: the automatically downloaded Ensembl annotation (or `--gtf`), registered for Stage 8.
- `haploblocks.<chrom>.tsv`: one row per haploblock.
- `cluster_memberships.<chrom>.tsv`: one row per represented sample haplotype and haploblock cluster. This is the independent cluster input for Stages 6 and 7.
- `sv_block_summary.<chrom>.tsv`: one row per overlapping VCF record and haploblock, including association counts but never duplicating a pair by passing cluster. This is the counting input for Stage 5.
- `sv_to_clusters.<chrom>.tsv`: one row per VCF record, haploblock, and cluster that passes the association threshold.

All paths are registered in `stage1_output/config.yaml`. Population labels remain independent of cluster inference: pass the Stage 0 table with `--sample-metadata` to normalize and publish it for Stages 4, 6, and 7. Useful method diagnostics are kept under `debug_and_qc/`; downloaded cluster files remain temporary and are removed after a successful run.

**Why we care:** Stage 1 creates the auditable allele-to-cluster mapping needed to ask whether the hash tags an SV, misses it within a cluster, or places the same resolved allele on several backgrounds.

The table keys and row meanings are part of the contract:

| Path key | One row per | Required identity columns |
|---|---|---|
| `sv_genotypes` | input VCF record | `sv_record_id` (unique); `sv_id`, `chrom`, `start`, `end`, and `sv_type` describe the record; sample GT columns follow the fixed metadata columns |
| `haploblocks` | haploblock | `haploblock_id`, `chrom`, `start`, `end` |
| `cluster_memberships` | complete sample haplotype assignment | `haploblock_id`, `sample_id`, `haplotype`, `cluster_id` |
| `sv_block_summary` | overlapping VCF-record–haploblock pair | `sv_record_id`, `haploblock_id`; this key is unique even when several clusters pass |
| `sv_to_clusters` | passing VCF-record–haploblock–cluster association | `sv_record_id`, `haploblock_id`, `cluster_id` |

The supplied `sample_metadata` table must provide `sample_id` and `population`. Stage 1 accepts original or canonical sample IDs, writes canonical IDs, removes metadata rows for samples absent from the VCF, and orders the result like the VCF. A `superpopulation` column may also be supplied, but downstream analyses must state explicitly which grouping they use.

### Downstream ownership and table grain

| Stage | Stage 1 inputs |
|---|---|
| 4 | `sv_genotypes` and `sample_metadata`; classify every `sv_record_id` once, before any haploblock join |
| 5 | `sv_block_summary` and `haploblocks`; count unique `sv_record_id, haploblock_id` pairs |
| 6 | `sv_genotypes`, `sv_block_summary`, `cluster_memberships`, and `sample_metadata`; test every record against every cluster in each overlapped block, including clusters that do not pass Stage 1's association threshold |
| 7 | `sv_genotypes`, `sv_block_summary`, `sv_to_clusters`, `cluster_memberships`, `sample_metadata`, and optional Stage 4 classification; measure carrier-cluster concentration and hash-group purity |
| 8 | Stage 6 `sv_cluster_summary`, Stage 7 `sv_hash_representation`, optional Stage 4 classifications, and the registered GTF; annotate one record–block candidate per Stage 6 summary row |

Use `sv_to_clusters` when the question is “which clusters pass the Stage 1 association rule for this SV?” It is intentionally thresholded and may contain several rows for a record. Use `sv_block_summary` with the distinct cluster IDs from `cluster_memberships` when the question is “which clusters could contain this SV?” That join preserves every overlapping record–block pair and every cluster defined in the block; it is the candidate universe used by Stages 6 and 7. Use `sv_genotypes` alone when every original VCF record must be retained independently of haploblock overlap or cluster assignment.

### Interpreting `sv_to_clusters.tsv`

Each row is an **SV–haploblock–cluster association**, not merely an overlap between an SV and an immutable block region. An SV can therefore have no passing cluster, one passing cluster, multiple passing clusters in the same block, or associations in more than one overlapped block. Absence from this table means that no evaluated cluster passed the configured probability and support thresholds; it does not prove that the SV is absent.

The main fields are:

| Field | Interpretation |
|---|---|
| `sv_record_id` | Unique, stable identifier for one input VCF record within the run, formatted `chrN_record_k`. Use it for joins and counts. |
| `sv_id` | Input VCF ID when it is short enough. It is a source label, not a unique key: callers can reuse it for records at the same coordinates. IDs longer than the configured limit are replaced with `SV_<chrom>_<start>_<end>_<type>_<hash>`. |
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

`pipeline/stage2_intersect.py` remains available for the separate descriptive boundary question. `position_class` uses exact overlap count: zero blocks is `outside_block`, one is `within_block`, and two or more is `boundary_crossing`. Proximity is retained separately as `near_boundary`, using `boundary_distance_bp`; a nearby SV is not mislabeled as crossing. It reuses Stage 1's `sv_genotypes` metadata tables, avoiding a second scan of the VCF, and reads Stage 1's generated config:

```bash
python claude-first-prototype/pipeline/stage2_intersect.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

This analysis is not required to establish cluster association. Because published blocks are contiguous regions, boundary results should be treated as a distinct spatial analysis rather than evidence that an SV belongs to a cluster.

**Why we care:** boundary-crossing events may reveal loci where the fixed block representation is structurally inappropriate, but this result does not by itself show that hashes improve SV prediction.

## Current Stage 4: population allele-frequency classification

`pipeline/stage4_classify_af.py` reads Stage 1's chromosome-specific `sv_genotypes` tables and normalized `sample_metadata.tsv`. It computes allele frequency for every population present in the metadata and classifies each SV once as `common`, `population_specific`, or `other`. It deliberately does not read haploblock cluster labels.

```bash
python claude-first-prototype/pipeline/stage4_classify_af.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

The stage writes `sv_af_classification.tsv` with one row per SV and population, `sv_classification.tsv` with one row per SV, and a `config.yaml` that carries the Stage 1 paths forward. Stage 6 can join the one-row-per-SV classification to `sv_block_summary` and independently compare it with `cluster_memberships`.

**Why we care:** population frequency is principally a control, but when joined to cluster prevalence it can distinguish ordinary ancestry structure from a population-enriched SV arising on an otherwise widespread haplotype background.

## Current Stage 5: per-haploblock SV-type enrichment

`pipeline/stage5_type_enrichment.py` reads Stage 1's chromosome-specific `sv_block_summary` and `haploblocks` tables. It counts each unique SV–haploblock pair once, builds the complete haploblock-by-SV-type grid, and estimates each cell's expected count from block length and the overall rate for that SV type. Two-sided Poisson p-values are corrected together with Benjamini–Hochberg FDR.

```bash
python claude-first-prototype/pipeline/stage5_type_enrichment.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

The output `stage5_output/sv_type_enrichment.tsv` contains observed and expected counts, p- and q-values, and a configurable significance flag for every block/type combination. The accompanying config carries the Stage 1 paths forward and registers this output as `paths.sv_type_enrichment`.

This is deliberately a readable first model with haploblock length as its only exposure adjustment. SNP density, callability, overdispersion, and minimum-count policies should be evaluated before treating significant cells as final biological results.

**Why we care:** a reproducible excess of one SV mechanism in a block could identify local sequence architecture that generates structural mutation, but it is optional because it does not directly test whether genomic hashes capture SV alleles.

## Rescoped Stage 6: population-conditioned SV–cluster association

`pipeline/stage6_cluster_association.py` evaluates each local cluster for every overlapping
SV–haploblock pair. It removes population means from cluster dosage and SV dosage before measuring
their association, then permutes SV genotypes within populations for an empirical null. This asks
whether the local hash contains information beyond population allele-frequency differences. The
outputs retain all tested cluster associations and a one-row-per-SV–block summary that distinguishes
cross-population-consistent tag candidates, population-dependent associations, exclusion signals,
and other detected associations.

The association table reports the number of called samples carrying each cluster, the SV carrier and
non-carrier counts among those samples, the carrier-rate difference, and whether the cluster is
carrier-enriched or carrier-depleted. Depleted clusters remain useful exclusion evidence but are not
reported as carrier tags.

**Why we care:** a cluster showing the same SV association across represented populations is a candidate portable tag, whereas an ancestry-only association is not.

The current implementation uses a per-SV–block maximum-statistic permutation and reports positive
carrier tags separately from negative exclusion signals. A small batched permutation run first
triages every pair; pairs passing that threshold receive an independent screening run, and promising
screening results receive a separate higher-resolution refinement before pair-level FDR correction.
Pairs whose observed effect is below `--min-abs-r` remain in the output but do not advance to the
costly screening or refinement stages. `permutations_used` records the size of the retained stage,
not the cumulative computation across discarded stages.
Population sign consistency remains an in-sample candidate screen rather than proof of portability;
sequencing-batch constraints and leave-one-population-out validation remain future refinements.
Genotypes are converted once per chromosome to a shared numeric dosage matrix, and Stage 6 reuses
pre-indexed block memberships and cluster matrices. Variants with the same call mask are analyzed
and permuted together using matrix operations rather than separate Python loops.

Run Stage 6 with Stage 4's carried-forward config when population classifications should remain
available to Stage 8:

```bash
python claude-first-prototype/pipeline/stage6_cluster_association.py \
  --config claude-first-prototype/stage4_output/config.yaml
```

## Rescoped Stage 7: hash representation audit and PCA QC

`pipeline/stage7_information_gain.py` measures how much local diplotype reduces uncertainty about
each SV and how often well-represented diplotypes contain both carriers and non-carriers. The first
quantity identifies SVs that are well tagged by existing hashes; the second identifies candidates
that add local information missing from those hashes. It also writes reusable PCA coordinate and
variance tables plus a population- and superpopulation-colored QC plot.

**Why we care:** Stage 7 is the direct “are haploblocks sufficient for SVs?” test: it measures whether each SV is concentrated in a small number of clusters and whether each cluster is internally homogeneous for SV carriage.

Raw in-sample information gain is not sufficient for that decision. The next revision should report
two directly interpretable quantities: carrier concentration across distinct clusters for each
resolved SV, and carrier/non-carrier purity among samples with each sufficiently represented cluster.
Because the VCF phase is not linked to the hash haplotype labels, a subdivision candidate is emitted
only when the same complete hash diplotype contains a configured minimum number of both carriers and
non-carriers. This is a count-supported candidate, not a replication result. Minimum
carrier/non-carrier counts, population- and batch-conditioned nulls, and resampling stability are
needed so rare calls or technical effects are not mistaken for hash failure. It should also flag
population-enriched SVs restricted within otherwise cosmopolitan clusters. PCA remains
batch/ancestry QC and is not evidence that the hash representation resolves SVs.

The proposed primary fields are `n_supported_carrier_clusters`, the fraction of supported carrier
evidence assigned to the top cluster, an effective carrier-cluster count derived from that evidence,
and per-cluster callable carrier/non-carrier counts and sample-level co-carriage purity. The
SV-centric concentration fields use Stage 1's allele-assignment evidence rather than counting both
haplotypes of every unphased carrier. Raw
information gain can remain a secondary descriptive column if its finite-sample limitations are
reported.

All-evidence concentration fields are diagnostic. Representation categories and population-context
flags use the separate standard-evidence top cluster, share, and effective cluster count so sparse
low-evidence assignments cannot create definitive multi-cluster or shared-background candidates.

```bash
python claude-first-prototype/pipeline/stage7_information_gain.py \
  --config claude-first-prototype/stage6_output/config.yaml \
  --threads 4
```

Stage 7 parses each chromosome's genotype table once and shares the resulting dosage matrix across
the representation audit and PCA. The diplotype summaries use block-level matrix operations, and
PCA uses an exact sample-by-sample eigendecomposition rather than decomposing the much wider
sample-by-SV matrix. `--threads N` processes independent chromosomes concurrently; choose a
modest value because each worker retains one chromosome dosage matrix in memory. For quick runs,
`--skip-information-gain` omits the secondary entropy calculation while preserving the direct
carrier-concentration and within-cluster subdivision results, and `--skip-pca` omits the structure
QC tables and plot entirely.

The primary outputs are `sv_carrier_cluster_summary.tsv`, `sv_cluster_purity.tsv`, and
`sv_hash_representation.tsv`. The existing `sv_haploblock_information.tsv` is retained as a secondary
description, and the PCA files remain QC outputs.

## Rescoped Stage 8: consequence-aware candidate annotation

`pipeline/stage8_candidate_annotation.py` joins Stage 6 candidates and Stage 7 representation patterns to the GTF registered by Stage 1 (or an optional `--gtf` override) and labels
consequences according to SV type. In particular, inversion breakpoint disruption is kept distinct
from genes merely contained in an inverted span. Association, representation, consequence, and
call-quality fields remain separate rather than being combined into a candidate score.

**Why we care:** Stage 8 adds gene, exon, and breakpoint context to the clearest examples of SVs captured or missed by the hash so those representation patterns can be interpreted at named loci.

Stage 8 should preserve the representation category—cross-population tag candidate, count-supported
hash-subdivision candidate, multi-cluster SV candidate, or population-enriched event on a shared
background—as a separate
field. Association, call quality, and functional annotation should remain separate evidence layers
rather than being collapsed into a single biological-importance score.

```bash
python claude-first-prototype/pipeline/stage8_candidate_annotation.py \
  --config claude-first-prototype/stage7_output/config.yaml
```

## Current testing status

Focused contract tests cover Stage 1's normalized downstream tables, metadata canonicalization, the distinction between exact boundary crossing and proximity, Stage 4 population classification, Stage 5's length-adjusted enrichment contract, population-conditioned Stage 6 controls, Stage 7 information gain/PCA tables, and SV-type-aware Stage 8 consequences. The earlier dbVar ingestion/QC tests described in the historical prototype no longer match the current pipeline.
