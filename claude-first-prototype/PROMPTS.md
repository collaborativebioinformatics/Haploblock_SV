# PROMPTS.md — Implementation handoff for the Haploblock_SV pipeline

Ready-to-use prompts for an agentic coding assistant, organized by pipeline stage and by day (see `README.md` for the full plan). Each prompt is self-contained enough to hand to a fresh agent session — it restates the relevant context rather than assuming prior conversation.

Shared context to paste into any session if it doesn't already have it: *"We're building a Python pipeline for a hackathon that studies how structural variants associate with haplotype-hash clusters inside haploblocks from data.haploblocks.org. The current input is a pre-merged 1000 Genomes ONT VCF containing DEL and INS calls. Stage 0 is a placeholder cohort-merging step being developed by Linh; the implemented Stage 1 reads the merged VCF directly and performs cluster-aware preprocessing. See README.md for the current contracts. Later prompts in this file preserve the original analysis ideas and are proposals, not implementation commitments."*

---

## Day 1 (Aug 25) — Stages 0–2

### Stage 0: Cohort SV merging (placeholder; Linh)

**Planned implementation:**
> Linh is working on Stage 0. It should download or accept a list of single-sample long-read SV VCFs, merge samples, and reconcile equivalent SV representations with `truvari collapse`. This is non-trivial because basecalling errors, mapping ambiguity around repeats, and caller differences can represent the same event differently. The eventual workflow will merge and collapse but will not run kanpig; it may therefore retain INV calls, while BND will probably remain unmerged.

**Current hackathon input:**
> Until Stage 0 exists, use `input/1kgp_ont_cohort.postfilter.full.vcf.gz`. It was generated previously with Sniffles, `bcftools merge`, `truvari collapse`, and kanpig regenotyping. Kanpig requires resolved sequences, so the present file includes DEL and INS but no INV or BND.

**Wire to next stage:**
> Stage 0's eventual primary output should be a merged cohort VCF accepted directly by Stage 1. It should also support replacement cohorts supplied as lists of single-sample VCFs.

### Stage 1: Cluster-aware preprocessing

**Implemented:**
> `pipeline/stage1_cluster_aware.py` reads the cohort VCF directly, downloads or accepts local per-block cluster membership files from data.haploblocks.org, preserves complete genotypes (with missing calls excluded from evidence), and infers SV-to-cluster associations using the other callable members of each sample's two clusters. It processes chromosomes independently and writes the downstream table `stage1_output/sv_to_clusters.tsv`, with useful diagnostics under `stage1_output/debug_and_qc/`.

**Interpretation:**
> A row links an SV to a haploblock **cluster**, not merely to the fixed block region. Associations require a cluster probability of at least 0.75 plus adaptive callable-haplotype support. Heterozygous assignments with posterior probability below 0.75 remain ambiguous. One SV may associate with multiple clusters in one block and may be evaluated independently in multiple overlapped blocks. See README.md for the field-by-field output contract.

**Wire to next stage:**
> Treat `sv_to_clusters.tsv` as the primary downstream input. Stage 2 is optional and answers the distinct spatial question of whether an SV is near or crosses a block boundary.

### Stage 2: Optional boundary classification

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

**Implement:**
> Write `pipeline/stage3_boundary_enrichment.py` that reads Stage 2's annotated SV table and tests whether `boundary_crossing` SVs are enriched relative to a null model. Implement the null as a circular rotation of SV positions within each chromosome (shift all SV coordinates on a chromosome by a random offset, wrapping around chromosome length, preserving intra-chromosome SV spacing/clustering while breaking their specific alignment to haploblock boundaries) — repeat for a configurable number of permutations (default 1000, seed from config.yaml) to build an empirical null distribution for the count of boundary-crossing SVs, and report an empirical p-value. Also implement a lightweight spatial-hotspot check: nearest-neighbor distance between consecutive SV breakpoints within each haploblock, compared to the expectation under a homogeneous Poisson process of the same intensity, reported as a z-score or KS-test p-value per block. Output a results table/report with both statistics.

**Test:**
> Write a pytest test for `stage3_boundary_enrichment.py` using two synthetic scenarios: (1) SVs placed exactly at random positions (should NOT show significant boundary enrichment, p > 0.05), and (2) SVs deliberately placed at every haploblock boundary (should show strong significant enrichment, p < 0.01). Fix the random seed and assert both outcomes.

**Wire to next stage:**
> No direct handoff needed for Stage 4 (it reads Stage 2's output independently), but confirm Stage 3's output is written to a path Stage 9 (integration) can later read for the summary report, and note that path in a shared `outputs.md` or in each stage's docstring.

### Stage 4: Common vs. population-specific SV classification

**Implemented:**
> `pipeline/stage4_classify_af.py` uses every SV in Stage 1's `sv_genotypes` tables and the population labels in its normalized `sample_metadata.tsv`. It classifies each SV once before any block or cluster join and writes both a per-population AF table and a one-row-per-SV classification table. The output classes are `common`, `population_specific`, and `other`, with `specific_to_population` and `other_reason` providing the relevant detail. The script never reads haploblock cluster labels, preserving the independence needed for Stage 6.

**Test:**
> Write a pytest test with a synthetic 3-SV, 2-population genotype matrix: one SV common to both populations, one private to population A only, one rare/absent everywhere. Assert the classification (`common`, `population_specific` with correct population, `other`) matches for all three at the default threshold.

**Wire to next stage:**
> Confirm Stage 6 (population-cluster correlation) and Stage 8 (DUP/INV overlay) both read the `sv_class`/`specific_to_population` columns added here, and that Stage 6 never reads data.haploblocks.org cluster labels from this script — only from its own separate input — to preserve the independence needed for a non-circular correlation test.

### Stage 5: Per-haploblock SV-type enrichment (with block-architecture offset)

**Implement:**
> Write `pipeline/stage5_type_enrichment.py` using Stage 1's `sv_block_summary` and `haploblocks` tables. Count each unique (`sv_id`, `haploblock_id`) pair once; do not count rows directly from `sv_to_clusters`, where one SV-block pair can have multiple passing clusters. Build a per-haploblock × SV-type count matrix. For each SV type, calculate the overall rate as the total number of SVs divided by total haploblock length, then calculate each block's expected count from its length. Use a Poisson test to compare observed and expected counts and apply Benjamini-Hochberg FDR correction across block × type tests. Output haploblock ID, SV type, observed count, expected count, p-value, FDR-adjusted q-value, and a q < 0.05 flag. This is a proposed analysis whose assumptions and power should be checked before implementation.

**Test:**
> Write a pytest test with synthetic data: several haploblocks of varying length with SV counts proportional to length (should show no significant enrichment after offset correction), plus one haploblock with an artificially inflated count for one SV type (should be flagged as significant after FDR correction). Assert the artificially-inflated block is flagged and the proportional ones are not.

**Wire to next stage:**
> Confirm Stage 9 (integration) reads this stage's flagged-haploblock table to build the "haploblocks prone to a particular SV type" section of the final report, and confirm the minimum-count exclusion fraction is surfaced there too as a documented limitation, not silently dropped.

---

## Day 3 (Aug 27) — Stages 6–8 (parallelizable)

### Stage 6: Population-cluster correlation

**Implement:**
> Write `pipeline/stage6_cluster_correlation.py` that reads Stage 4's classified SV table and computes, per haploblock, a population-specific-SV density (e.g. count of population-specific SVs / total SVs, or / block length). Use the cluster memberships downloaded and normalized by Stage 1 to derive the separate SNV-haplotype cluster quantity to compare against. Compute a Spearman correlation between the two per-block quantities, with a permutation-based p-value (shuffle block labels, default 1000 permutations). Output the correlation coefficient, p-value, and a per-block table for plotting. The exact cluster differentiation summary still needs biological definition.

**Test:**
> Write a pytest test with synthetic data where population-specific-SV density is constructed to be a noisy linear function of a synthetic cluster-differentiation score; assert the computed Spearman correlation is positive and significant (p < 0.05) at a fixed seed, and a shuffled-label negative control on the same data is not significant.

**Wire to next stage:**
> Confirm Stage 9 reads this stage's correlation coefficient/p-value and per-block table to render the population-correlation section of the report; note in a comment whether the sign/magnitude found matches or contradicts H3's expectation, so the report can state the finding plainly.

### Stage 7: SV-based population structure reconstruction

**Implement:**
> Write `pipeline/stage7_sv_clustering.py` that builds a per-sample × per-haploblock SV presence/absence (or dosage) matrix, runs PCA and UMAP with fixed seeds, clusters the reduced embedding, and compares the assignments with (a) the populations supplied by `sample_metadata.tsv` and (b) data.haploblocks.org's hash-based clusters. Output PNG PCA/UMAP plots and cluster assignments; additional agreement statistics remain optional pending a precise definition of the cross-block cluster labels.

**Test:**
> Write a pytest test with a small synthetic SV matrix constructed so that two groups of samples have clearly distinct SV profiles (e.g. disjoint sets of population-specific SVs); assert the resulting clustering recovers the two groups with ARI > 0.8 against the known synthetic group labels, at a fixed seed.

**Wire to next stage:**
> Confirm Stage 9 can read the embedding and cluster assignments to render scatter plots colored by the populations from `sample_metadata.tsv` and, separately, by an appropriately defined haploblocks.org cluster label. If agreement statistics are retained, report differences plainly.

### Stage 8: Duplication/inversion gene overlay

**Implement:**
> For a future cohort VCF containing DUP and/or INV records, write `pipeline/stage8_dup_inv_overlay.py` to identify recurrent and population-specific calls and overlap them against a gene annotation GTF. Output a ranked table of stand-out DUP/INV calls with overlapping gene names. The current kanpig-regenotyped input contains only DEL and INS, so this stage cannot run meaningfully on it. A selection-scan overlay is deferred to possible future work rather than included in this stage.

**Test:**
> Write a pytest test with a synthetic DUP/INV table and a tiny synthetic GTF (2-3 genes at known coordinates); assert the gene-overlap output correctly attributes each synthetic SV to the gene(s) it overlaps and correctly excludes the ones that don't overlap any gene.

**Wire to next stage:**
> Confirm Stage 9 reads this stage's ranked DUP/INV table for the report section specifically called out for Maria (duplications) and Alistair (inversions), and that the table is sorted/filterable by SV type so each of them can look at just their type of interest.

---

## Final section: integration, documentation, and demo

**Integration test (future, after stages are selected and implemented):**
> Write an end-to-end integration test on a single chromosome that runs only the stages retained in the final analysis plan, asserts each expected output is non-empty, and validates the final summary contract. Stage 0 should be included once Linh's merging workflow is implemented.

**Documentation:**
> Review `README.md`'s pipeline overview table against the actual `pipeline/` scripts as implemented, and update either the README or the code so stage names, inputs, and outputs match exactly. Add a `--help` docstring/argparse description to every stage script consistent with the README's one-line purpose for that stage.

**Demo / report:**
> If the proposed downstream stages are retained, write `pipeline/stage9_integrate.py` to aggregate their outputs into a per-haploblock summary and render selected figures in a notebook or static HTML report. A web application is outside the current scope and can be reconsidered as future work.
