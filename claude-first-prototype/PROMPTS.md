# PROMPTS.md — Implementation handoff for the Haploblock_SV pipeline

Ready-to-use prompts for an agentic coding assistant, organized by pipeline stage and by day (see `README.md` for the full plan). Each prompt is self-contained enough to hand to a fresh agent session — it restates the relevant context rather than assuming prior conversation.

Shared context to paste into any session if it doesn't already have it: *"We're building a Python pipeline for a hackathon that studies how structural variants associate with haplotype-hash clusters inside haploblocks from data.haploblocks.org. The current input is a pre-merged 1000 Genomes ONT VCF containing DEL and INS calls. Stage 0 is a placeholder cohort-merging step being developed by Linh; the implemented Stage 1 reads the merged VCF directly and performs cluster-aware preprocessing. See README.md for the current contracts. Later prompts in this file preserve the original analysis ideas and are proposals, not implementation commitments."*

Shared biological decision to paste into downstream sessions: *"Treat this as an evaluation of the
haploblock representation. Ask whether carriers of one resolved SV are concentrated in a small number
of genomic-hash clusters and whether each cluster is homogeneous for SV carriage. Population structure
is a control and biological context, not an alternative SV-prediction model or the endpoint."*

---

## Day 1 (Aug 25) — Stages 0–2

### Stage 0: Cohort SV merging (placeholder; Linh)

**Why we care:**
> Equivalent biological SV alleles must be represented consistently before multi-cluster carriage can be interpreted as recurrence or haplotype history rather than call fragmentation.

**Planned implementation:**
> Linh is working on Stage 0. It should download or accept a list of single-sample long-read SV VCFs, merge samples, and reconcile equivalent SV representations with `truvari collapse`. This is non-trivial because basecalling errors, mapping ambiguity around repeats, and caller differences can represent the same event differently. The eventual workflow will merge and collapse but will not run kanpig; it may therefore retain INV calls, while BND will probably remain unmerged.

**Current hackathon input:**
> Until Stage 0 exists, use `input/1kgp_ont_cohort.postfilter.full.vcf.gz`. It was generated previously with Sniffles, `bcftools merge`, `truvari collapse`, and kanpig regenotyping. Kanpig requires resolved sequences, so the present file includes DEL and INS but no INV or BND.

**Wire to next stage:**
> Stage 0's eventual primary output should be a merged cohort VCF accepted directly by Stage 1. It should also support replacement cohorts supplied as lists of single-sample VCFs.

### Stage 1: Cluster-aware preprocessing

**Why we care:**
> This stage creates the auditable allele-to-cluster mapping needed to test what the genomic hash captures, misses, or assigns to several backgrounds.

**Implemented:**
> `pipeline/stage1_cluster_aware.py` reads the cohort VCF directly, downloads or accepts local per-block cluster membership files from data.haploblocks.org, preserves complete genotypes (with missing calls excluded from evidence), and infers SV-to-cluster associations using the other callable members of each sample's two clusters. It processes chromosomes independently and writes the downstream table `stage1_output/sv_to_clusters.tsv`, with useful diagnostics under `stage1_output/debug_and_qc/`.

**Interpretation:**
> A row links an SV to a haploblock **cluster**, not merely to the fixed block region. Associations require a cluster probability of at least 0.75 plus adaptive callable-haplotype support. Heterozygous assignments with posterior probability below 0.75 remain ambiguous. One SV may associate with multiple clusters in one block and may be evaluated independently in multiple overlapped blocks. See README.md for the field-by-field output contract.

**Wire to next stage:**
> Treat `sv_to_clusters.tsv` as the primary downstream input. Stage 2 is optional and answers the distinct spatial question of whether an SV is near or crosses a block boundary.

### Stage 2: Optional boundary classification

**Why we care:**
> Boundary-crossing SVs can expose loci where a fixed block is a poor representation, but boundary position alone does not establish that haploblock clusters are useful.

**Implemented:**
> `pipeline/stage2_intersect.py` may be run from Stage 1's generated config. Its `position_class` is based on exact overlap count (`outside_block`, `within_block`, or `boundary_crossing`); boundary proximity is a separate `near_boundary` boolean. Keep this descriptive result separate from cluster association: overlap with an immutable block region does not show that an SV belongs to a haplotype cluster.

**Use:**
> Retain Stage 2 only if boundary-crossing SVs remain a biological question of interest. It is not required for cluster-aware downstream analyses.

**Wire to next stage:**
> Only the proposed boundary-enrichment analysis in Stage 3 depends on this output. Other downstream work should start from Stage 1's cluster-aware table and, when required, the original genotype information.

---

## Day 2 (Aug 26) — Stages 3–5

The remaining prompts are preserved from the original hackathon plan. They are candidate analyses that need biological review and adaptation to the cluster-aware Stage 1 output before implementation.

### Stage 3: Boundary enrichment test

**Why we care:**
> A reproducible boundary excess could suggest that recombination-defined block boundaries coincide with structural instability; otherwise this stage should remain secondary to the cluster-representation questions.

**Implement:**
> Write `pipeline/stage3_boundary_enrichment.py` that reads Stage 2's annotated SV table and tests whether `boundary_crossing` SVs are enriched relative to a null model. Implement the null as a circular rotation of SV positions within each chromosome (shift all SV coordinates on a chromosome by a random offset, wrapping around chromosome length, preserving intra-chromosome SV spacing/clustering while breaking their specific alignment to haploblock boundaries) — repeat for a configurable number of permutations (default 1000, seed from config.yaml) to build an empirical null distribution for the count of boundary-crossing SVs, and report an empirical p-value. Also implement a lightweight spatial-hotspot check: nearest-neighbor distance between consecutive SV breakpoints within each haploblock, compared to the expectation under a homogeneous Poisson process of the same intensity, reported as a z-score or KS-test p-value per block. Output a results table/report with both statistics.

**Test:**
> Write a pytest test for `stage3_boundary_enrichment.py` using two synthetic scenarios: (1) SVs placed exactly at random positions (should NOT show significant boundary enrichment, p > 0.05), and (2) SVs deliberately placed at every haploblock boundary (should show strong significant enrichment, p < 0.01). Fix the random seed and assert both outcomes.

**Wire to next stage:**
> No direct handoff needed for Stage 4 (it reads Stage 2's output independently), but confirm Stage 3's output is written to a path Stage 9 (integration) can later read for the summary report, and note that path in a shared `outputs.md` or in each stage's docstring.

### Stage 4: Common vs. population-specific SV classification

**Why we care:**
> Population frequency controls ancestry confounding and, when joined to cluster prevalence, can reveal a population-enriched SV occurring within an otherwise cosmopolitan haplotype background.

**Implemented:**
> `pipeline/stage4_classify_af.py` uses every SV in Stage 1's `sv_genotypes` tables and the population labels in its normalized `sample_metadata.tsv`. It classifies each SV once before any block or cluster join and writes both a per-population AF table and a one-row-per-SV classification table. The output classes are `common`, `population_specific`, and `other`, with `specific_to_population` and `other_reason` providing the relevant detail. The script never reads haploblock cluster labels, preserving the independence needed for Stage 6.

**Test:**
> Write a pytest test with a synthetic 3-SV, 2-population genotype matrix: one SV common to both populations, one private to population A only, one rare/absent everywhere. Assert the classification (`common`, `population_specific` with correct population, `other`) matches for all three at the default threshold.

**Wire to next stage:**
> Confirm Stage 6 (population-cluster correlation) and Stage 8 (DUP/INV overlay) both read the `sv_class`/`specific_to_population` columns added here, and that Stage 6 never reads data.haploblocks.org cluster labels from this script — only from its own separate input — to preserve the independence needed for a non-circular correlation test.

### Stage 5: Per-haploblock SV-type enrichment (with block-architecture offset)

**Why we care:**
> Type enrichment can nominate block-local sequence architecture that generates SVs, but it is optional because it does not directly test whether genomic hashes capture individual SV alleles.

**Implemented:**
> `pipeline/stage5_type_enrichment.py` uses Stage 1's `sv_block_summary` and `haploblocks` tables. It counts each unique (`sv_id`, `haploblock_id`) pair once, builds the complete haploblock × SV-type matrix, calculates length-adjusted expected counts, and applies two-sided Poisson tests with Benjamini-Hochberg correction across all cells. It writes observed and expected counts, p- and q-values, and the configurable significance flag to `sv_type_enrichment.tsv`, then registers that path in a carried-forward config. The length-only Poisson model is the current prototype; callability, SNP density, and overdispersion remain candidates for model refinement.

**Test:**
> Write a pytest test with synthetic data: several haploblocks of varying length with SV counts proportional to length (should show no significant enrichment after offset correction), plus one haploblock with an artificially inflated count for one SV type (should be flagged as significant after FDR correction). Assert the artificially-inflated block is flagged and the proportional ones are not.

**Wire to next stage:**
> Stage 9 (integration) should read `paths.sv_type_enrichment` from the carried-forward config and use its flagged cells for the "haploblocks prone to a particular SV type" section. Report that this prototype adjusts for block length only; do not imply that SNP density, callability, or overdispersion have already been modeled.

---

## Day 3 (Aug 27) — Rescoped Stages 6–8 (parallelizable)

The original proposals below were refined after biological review. Population structure is an
important control and descriptive result, but SV-based ancestry plots and population-label
agreement are already well established. The revised analyses therefore ask a more specific
question: what information about SV carriage is captured, or missed, by the local SNV-derived
haploblock clusters? This keeps the useful population comparisons while making the main outputs
interpretable at individual loci. Functional annotation is then applied to the most informative
SVs instead of treating any gene overlap as equivalent evidence.

### Stage 6: Population-conditioned SV–cluster association

**Why we care:**
> A hash cluster is useful when it acts as a portable proxy for an SV allele in samples or populations without long-read genotypes; association caused only by shared ancestry or sequencing batch does not provide that utility.

**Implement:**
> Write `pipeline/stage6_cluster_association.py` using Stage 1's `sv_genotypes`,
> `sv_block_summary`, `cluster_memberships`, and `sample_metadata` tables. For every overlapping
> SV–haploblock pair, test whether dosage of each local SNV-derived cluster predicts SV dosage
> after population means are removed. Estimate an empirical p-value by permuting SV genotypes
> within populations, thereby preserving population-specific allele frequencies. Report carrier
> rates with and without the cluster, the population-adjusted correlation, FDR, the number of
> individually informative populations, and directional consistency across them. Classify the
> strongest association per SV–block as a cross-population-consistent tag candidate, population-dependent association,
> cluster association, or no detected cluster signal. Stage 6 should read normalized cluster tables
> from Stage 1's config; cluster fetching and normalization remain Stage 1 responsibilities rather
> than being duplicated in Stage 0 or Stage 6.

**Required refinement before biological claims:**
> Treat selection of the best cluster as part of the test by permuting the maximum absolute effect
> across all eligible clusters for each SV–block pair. Constrain permutations jointly by population
> and sequencing batch when batch is available, use enough adaptive permutations to resolve the
> desired FDR, and report positive carrier enrichment separately from negative cluster exclusion.
> For candidate portable tags, emit per-population carrier rates/effects and leave-one-population-out
> validation rather than defining portability from sign consistency in the discovery samples alone.
> Use a modest maximum-statistic screening run for every SV–block pair, then discard its p-value and
> run an independent high-resolution permutation set for promising pairs before pair-level FDR.

**Test:**
> Use synthetic populations in which one cluster predicts an SV within both populations, plus an
> SV whose apparent association is caused only by different population frequencies. Assert that
> the within-population permutation identifies the portable local association and does not promote
> the ancestry-only signal.

**Wire to next stage:**
> Stage 8 should use the per-SV association and pattern fields to prioritize functional annotation.
> Stage 9 should summarize counts and examples of each association pattern rather than elevate one
> genome-wide correlation as the biological result.

### Stage 7: Hash representation audit and population-structure QC

**Why we care:**
> This stage decides whether haploblocks sufficiently resolve SV carriage: the same SV appearing across many standard-evidence clusters and a count-supported carrier/non-carrier split within one hash diplotype are the two primary ways the hash can miss SV information.

**Implement:**
> Write `pipeline/stage7_information_gain.py` that measures, per SV–haploblock pair, how much local
> haploblock diplotype reduces uncertainty about SV carriage and how often sufficiently represented
> diplotypes still contain both carriers and non-carriers. Aggregate these quantities per block to
> identify blocks where SVs are well tagged and blocks where SVs add information missing from the
> SNV-derived hashes. Also produce a genome-wide SV PCA as a QC view, saving coordinate and
> explained-variance tables as well as a PNG colored by `population` and, when available,
> `superpopulation`. PCA is descriptive QC rather than the main biological result; UMAP and
> unsupervised population clustering are omitted until they answer a defined downstream question.

**Required refinement before biological claims:**
> Replace raw in-sample normalized information gain as the headline result with two direct summaries:
> carrier concentration across distinct clusters for each resolved SV, and carrier/non-carrier purity
> among samples with each sufficiently represented cluster. Because VCF phase is not linked to hash
> haplotype labels, require the same complete diplotype to contain enough carriers and non-carriers
> before emitting a hash-subdivision candidate. Require minimum carrier, non-carrier, and cluster
> counts, assess stability by resampling, and use population- and batch-conditioned nulls. The primary
> output should state where haploblocks are sufficient SV tags and where direct SV genotypes add
> resolution, rather than emphasizing an abstract information score.

**Primary output fields:**
> For each resolved SV–block pair, report the number of supported carrier clusters, the fraction of
> supported carrier evidence assigned to the top cluster, and an effective carrier-cluster count.
> For each tested cluster, report callable carrier and non-carrier counts and sample-level co-carriage
> purity; separately report the number of mixed complete diplotypes meeting the configured carrier
> and non-carrier count threshold. Treat these as subdivision candidates rather than replicated results. Use Stage 1's
> allele-assignment evidence for SV-centric cluster concentration rather than naively counting both
> haplotypes of an unphased carrier. Retain normalized information gain only as a
> secondary descriptive field with its finite-sample limitation documented.
> Report all-evidence breadth as a diagnostic, but derive representation and population-context
> categories only from standard-evidence carrier assignments.

**Add explicit representation-pattern outputs:**
> After multiallelic/VNTR alleles have been resolved at the locus level, flag (1) the same allele on
> multiple divergent carrier clusters and (2) a Stage 4 population-enriched SV whose carrier cluster
> is cosmopolitan but whose SV carriage is population restricted within that cluster. Require
> replication across batches or read/sequence review before interpreting the first pattern as
> recurrence, gene conversion, recombination, or allele age.

**Test:**
> Construct one SV perfectly tagged by local diplotype and one that segregates within every
> sufficiently represented diplotype. Assert high information gain for the first and a high
> within-diplotype mixed fraction for the second. Confirm that PCA tables are deterministic.

**Wire to next stage:**
> Stage 9 should report blocks and SVs at both extremes: portable/taggable variation and variation
> that subdivides existing hashes. Population-colored PCA remains a QC figure.
> Run Stage 7 from Stage 6's carried-forward config so `sv_hash_representation` and the Stage 6
> association paths are both available to Stage 8.

### Stage 8: Consequence-aware candidate annotation

**Why we care:**
> Gene, exon, and breakpoint annotation makes the clearest examples of information captured or missed by the hash interpretable at specific loci without changing the core representation question.

**Implement:**
> Write `pipeline/stage8_candidate_annotation.py` to annotate Stage 6 candidates against a GTF.
> Interpret overlap by SV type: exon and transcript loss for DEL, complete or partial gene overlap
> for DUP, breakpoint disruption separately from genes merely contained inside INV spans, and
> breakpoint context for INS. Rank candidates using transparent evidence fields including Stage 6
> association pattern/strength, population classification when available, call precision/filter,
> and consequence class. Run on DEL/INS now; preserve the same contract for future DUP/INV inputs.
> Constraint, regulatory, repeat-context, and well-controlled selection evidence
> can be joined later without reducing them to a single proximity boolean.

**Required refinement before biological claims:**
> Carry the representation pattern from Stages 6–7 into the candidate table: cross-population tag
> candidate, count-supported hash-subdivision candidate, standard-evidence multi-cluster SV
> candidate, or population-enriched event on a shared background. Keep call quality, association,
> transcript consequence, and
> available annotation evidence as separate sortable fields. Do not let a single additive score hide
> whether a candidate is statistically supported, technically credible, or merely gene-overlapping.

**Test:**
> Use a tiny synthetic GTF with exons and genes around DEL, INS, DUP, and INV examples. Assert that
> exon loss, insertion breakpoint context, duplication overlap, and inversion breakpoint disruption
> are distinguished, and that a gene only contained within an inversion is not described as a
> breakpoint-disrupted gene.

**Wire to next stage:**
> Stage 9 should expose a filterable candidate table by SV type and evidence category, with DEL/INS
> available for the present callset and DUP/INV views enabled when those calls are available.
> Run Stage 8 from Stage 7's carried-forward config so representation patterns are joined to the
> Stage 6 candidates.

---

## Final section: integration, documentation, and demo

**Integration test (future, after stages are selected and implemented):**
> Write an end-to-end integration test on a single chromosome that runs only the stages retained in the final analysis plan, asserts each expected output is non-empty, and validates the final summary contract. Stage 0 should be included once Linh's merging workflow is implemented.

**Documentation:**
> Review `README.md`'s pipeline overview table against the actual `pipeline/` scripts as implemented, and update either the README or the code so stage names, inputs, and outputs match exactly. Add a `--help` docstring/argparse description to every stage script consistent with the README's one-line purpose for that stage.

**Demo / report:**
> If the proposed downstream stages are retained, write `pipeline/stage9_integrate.py` to aggregate their outputs into a per-haploblock summary and render selected figures in a notebook or static HTML report. Organize the report around where SVs are concentrated in a few clusters, where the same SV spans many clusters, and where individual clusters split into carriers and non-carriers rather than around one panel per metric. This matters because the integration stage must state where haploblocks do and do not resolve SV carriage, not merely demonstrate that every script ran. A web application is outside the current scope and can be reconsidered as future work.
