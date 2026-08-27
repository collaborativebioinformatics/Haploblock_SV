# PROMPTS.md — Implementation handoff for the Haploblock_SV pipeline

Ready-to-use prompts for an agentic coding assistant, organized by pipeline stage and by day (see `README.md` for the full plan). Each prompt is self-contained enough to hand to a fresh agent session — it restates the relevant context rather than assuming prior conversation.

Shared context to paste into any session if it doesn't already have it: *"We're building a Python pipeline for a hackathon that studies structural variants (SVs: DEL/DUP/INV/INS) within 'haploblocks' — LD-defined haplotype-hash regions from data.haploblocks.org. SV calls come from dbVar study nstd152 (Chaisson et al. 2019, 1000 Genomes haplotype-resolved SVs). The pipeline is a sequence of independently-runnable, parameterized stages under `pipeline/`, each a script or module taking input paths and config values as CLI args or function args — never hardcoded paths. See README.md for the full stage list and PROMPTS.md for the stage you're implementing."*

This is an example change for Maria to see.

---

## Day 1 (Aug 25) — Stages 0–2

### Stage 0: Data ingestion & harmonization

**Implement:**
> Write `pipeline/stage0_ingest.py`, a standalone script that: (1) downloads or accepts a local path to the dbVar nstd152 VCF (structural variants), (2) downloads or accepts a local path to haploblock BED/metadata and population-cluster labels from data.haploblocks.org, (3) downloads or accepts a local path to 1000 Genomes sample→superpopulation metadata (AFR/AMR/EAS/EUR/SAS), (4) confirms all inputs are on the same genome build (liftover to GRCh38 with `pyliftover` if not, logging a warning), and (5) writes one shared `config.yaml` capturing genome build, AF thresholds to be used later, boundary-distance threshold N (bp), random seeds for permutation/UMAP, and file paths — this config is read by every later stage, so its schema is the contract for the rest of the pipeline. Take all input paths/URLs as CLI args with sensible defaults, not hardcoded. If any real data source is slow or its schema is unclear, fall back to generating a small synthetic dataset (a few hundred fake SVs across a few dozen fake haploblocks with plausible fields) so downstream stages are never blocked — log clearly when running on synthetic vs. real data.

**Test:**
> Write a small pytest test for `stage0_ingest.py` that runs it in synthetic-data mode (no network access) and asserts: `config.yaml` is created and parses as valid YAML with all expected keys (genome build, thresholds, seeds, paths), the SV table and haploblock table are non-empty pandas DataFrames with the expected columns, and sample metadata has superpopulation labels drawn from the standard five values.

**Wire to next stage:**
> Confirm `stage0_ingest.py`'s output file paths (SV table, haploblock table, sample metadata, config.yaml) match exactly what `pipeline/stage1_qc.py` expects to read — if Stage 1 doesn't exist yet, write down the expected input contract (file names, formats, columns) as a short comment block at the top of `stage0_ingest.py` so Stage 1 can be implemented against it independently.

### Stage 1: QC & normalization

**Implement:**
> Write `pipeline/stage1_qc.py`, a standalone script taking the raw SV table and haploblock table (paths from Stage 0's config.yaml) and: filtering SVs to PASS/high-confidence calls only, dropping SVs below/above configurable size thresholds (read from config.yaml), normalizing breakpoint coordinates (left-aligned, 0-based BED-style), deduplicating exact-duplicate calls, and validating the haploblock BED is sorted and non-overlapping (raise a clear error listing offending rows if not, rather than silently continuing). Output a cleaned SV table and cleaned haploblock table as parquet or TSV, plus a short QC report (counts before/after each filter) as a JSON or text file.

**Test:**
> Write a pytest test for `stage1_qc.py` using a small hand-crafted SV table containing: one PASS call, one low-confidence call that should be dropped, one duplicate call, and one call below the minimum size threshold. Assert the output table retains only the one valid call and the QC report's counts match expectations exactly.

**Wire to next stage:**
> Update `pipeline/stage2_intersect.py` (or its docstring/contract if not yet written) to read Stage 1's cleaned SV table and cleaned haploblock table as its inputs, and confirm the coordinate convention (0-based BED-style, as fixed in Stage 1) is documented and consistent between the two scripts.

### Stage 2: SV × haploblock intersection

**Implement:**
> Write `pipeline/stage2_intersect.py` using `pybedtools` to intersect the QC'd SV table (as BED-like intervals) against the QC'd haploblock BED. For each SV, add a `position_class` column: `within_block` (fully contained), `boundary_crossing` (either breakpoint within N bp of a block edge — N read from config.yaml — but not fully outside), or `outside_block`. Also add a `haploblock_id` column (or a list, if an SV can touch multiple blocks) linking each SV back to the block(s) it overlaps. Output one annotated SV table (parquet/TSV) that is the shared input for all downstream stages. Take all paths as CLI args, not hardcoded.

**Test:**
> Write a pytest test for `stage2_intersect.py` with a synthetic haploblock (e.g. chr1:1000-2000) and three synthetic SVs: one fully inside (expect `within_block`), one straddling the 2000 boundary within N bp (expect `boundary_crossing`), and one far away on chr1 (expect `outside_block`). Assert the output `position_class` column matches exactly for all three.

**Wire to next stage:**
> Confirm the annotated SV table's schema (columns: SV id, type, chrom, start, end, haploblock_id, position_class, sample genotype columns) is exactly what Stage 3 (boundary enrichment) and Stage 4 (AF classification) both need, since they read this same table independently — document the schema as a comment block at the top of `stage2_intersect.py`.

---

## Day 2 (Aug 26) — Stages 3–5

### Stage 3: Boundary enrichment test

**Implement:**
> Write `pipeline/stage3_boundary_enrichment.py` that reads Stage 2's annotated SV table and tests whether `boundary_crossing` SVs are enriched relative to a null model. Implement the null as a circular rotation of SV positions within each chromosome (shift all SV coordinates on a chromosome by a random offset, wrapping around chromosome length, preserving intra-chromosome SV spacing/clustering while breaking their specific alignment to haploblock boundaries) — repeat for a configurable number of permutations (default 1000, seed from config.yaml) to build an empirical null distribution for the count of boundary-crossing SVs, and report an empirical p-value. Also implement a lightweight spatial-hotspot check: nearest-neighbor distance between consecutive SV breakpoints within each haploblock, compared to the expectation under a homogeneous Poisson process of the same intensity, reported as a z-score or KS-test p-value per block. Output a results table/report with both statistics.

**Test:**
> Write a pytest test for `stage3_boundary_enrichment.py` using two synthetic scenarios: (1) SVs placed exactly at random positions (should NOT show significant boundary enrichment, p > 0.05), and (2) SVs deliberately placed at every haploblock boundary (should show strong significant enrichment, p < 0.01). Fix the random seed and assert both outcomes.

**Wire to next stage:**
> No direct handoff needed for Stage 4 (it reads Stage 2's output independently), but confirm Stage 3's output is written to a path Stage 9 (integration) can later read for the summary report, and note that path in a shared `outputs.md` or in each stage's docstring.

### Stage 4: Common vs. population-specific SV classification

**Implement:**
> Write `pipeline/stage4_classify_af.py` that reads Stage 2's annotated SV table (with per-sample genotype columns) and Stage 0's sample→superpopulation metadata, computes allele frequency per SV within each of the five standard 1000 Genomes superpopulations (AFR/AMR/EAS/EUR/SAS — NOT the data.haploblocks.org clusters, to avoid circularity with Stage 6), and classifies each SV as `common` (AF above a configurable threshold, default 0.05, in ≥2 superpopulations) or `population_specific` (AF above threshold in exactly one superpopulation, near-zero elsewhere) or `other` (doesn't meet either definition — log the fraction falling here). Add `sv_class` and `specific_to_population` (nullable) columns to the SV table and write it out.

**Test:**
> Write a pytest test with a synthetic 3-SV, 2-superpopulation genotype matrix: one SV common to both populations, one private to population A only, one rare/absent everywhere. Assert the classification (`common`, `population_specific` with correct population, `other`) matches for all three at the default threshold.

**Wire to next stage:**
> Confirm Stage 6 (population-cluster correlation) and Stage 8 (DUP/INV overlay) both read the `sv_class`/`specific_to_population` columns added here, and that Stage 6 never reads data.haploblocks.org cluster labels from this script — only from its own separate input — to preserve the independence needed for a non-circular correlation test.

### Stage 5: Per-haploblock SV-type enrichment (with block-architecture offset)

**Implement:**
> Write `pipeline/stage5_type_enrichment.py` that builds a per-haploblock × SV-type count matrix from Stage 2's output. For each SV type, calculate the overall SV rate across all haploblocks as the total number of SVs of that type divided by the total haploblock length. For each haploblock, calculate the expected number of SVs of that type based on its length and the overall SV rate. Use a Poisson test to determine whether the observed SV count significantly deviates from the expected count. Apply Benjamini-Hochberg FDR correction across all haploblock × SV-type tests. Flag combinations with an FDR-adjusted q-value below 0.05.


**Test:**
> Write a pytest test with synthetic data: several haploblocks of varying length with SV counts proportional to length (should show no significant enrichment after offset correction), plus one haploblock with an artificially inflated count for one SV type (should be flagged as significant after FDR correction). Assert the artificially-inflated block is flagged and the proportional ones are not.

**Wire to next stage:**
> Confirm Stage 9 (integration) reads this stage's flagged-haploblock table to build the "haploblocks prone to a particular SV type" section of the final report, and confirm the minimum-count exclusion fraction is surfaced there too as a documented limitation, not silently dropped.

---

## Day 3 (Aug 27) — Stages 6–8 (parallelizable)

### Stage 6: Population-cluster correlation

**Implement:**
> Write `pipeline/stage6_cluster_correlation.py` that reads Stage 4's classified SV table and computes, per haploblock, a population-specific-SV density (e.g. count of population-specific SVs / total SVs, or / block length). Separately read data.haploblocks.org's per-haploblock cluster labels/differentiation metric (from Stage 0's output) and compute a correlation (Spearman, since the relationship needn't be linear) between the two per-block quantities, with a permutation-based p-value (shuffle block labels, default 1000 permutations). Output the correlation coefficient, p-value, and a per-block table for plotting.

**Test:**
> Write a pytest test with synthetic data where population-specific-SV density is constructed to be a noisy linear function of a synthetic cluster-differentiation score; assert the computed Spearman correlation is positive and significant (p < 0.05) at a fixed seed, and a shuffled-label negative control on the same data is not significant.

**Wire to next stage:**
> Confirm Stage 9 reads this stage's correlation coefficient/p-value and per-block table to render the population-correlation section of the report; note in a comment whether the sign/magnitude found matches or contradicts H3's expectation, so the report can state the finding plainly.

### Stage 7: SV-based population structure reconstruction

**Implement:**
> Write `pipeline/stage7_sv_clustering.py` that builds a per-sample × per-haploblock SV presence/absence (or dosage) matrix from Stage 2's output and genotype columns, runs PCA (via scikit-learn) and UMAP (via umap-learn, fixed `random_state` from config.yaml) for dimensionality reduction, clusters the reduced embedding (e.g. k-means or HDBSCAN), and computes Adjusted Rand Index between the resulting clusters and (a) 1000 Genomes superpopulation labels and (b) data.haploblocks.org's own hash-based clusters. Output the embedding coordinates, cluster assignments, and both ARI scores.

**Test:**
> Write a pytest test with a small synthetic SV matrix constructed so that two groups of samples have clearly distinct SV profiles (e.g. disjoint sets of population-specific SVs); assert the resulting clustering recovers the two groups with ARI > 0.8 against the known synthetic group labels, at a fixed seed.

**Wire to next stage:**
> Confirm Stage 9 reads both ARI scores and the embedding coordinates to render a scatter plot (colored by superpopulation, and separately by haploblocks.org cluster) in the final report, and flag prominently if the ARI against superpopulations is much higher than against haploblocks.org clusters (or vice versa) since that's a key H6 finding either way.

### Stage 8: Duplication/inversion gene & selection overlay

**Implement:**
> Write `pipeline/stage8_dup_inv_overlay.py` that filters Stage 4's classified SV table to DUP and INV calls, identifies recurrent ones (present in multiple samples) and population-specific ones, and overlaps them against a gene annotation GTF (via `pyranges`/`gtfparse`) to report which genes each recurrent/population-specific DUP or INV overlaps. Output a ranked table (by recurrence count, then by population specificity) of "stand-out" DUP/INV calls with overlapping gene names — this is the core, always-in-scope part. As a clearly separate optional function `add_selection_scan_overlay()` (Stage 8b, stretch goal, skip if out of time), cross-reference the same calls against a selection-scan region file (e.g. PopHumanScan BED) if one is available, adding a `near_selection_signal` boolean column.

**Test:**
> Write a pytest test with a synthetic DUP/INV table and a tiny synthetic GTF (2-3 genes at known coordinates); assert the gene-overlap output correctly attributes each synthetic SV to the gene(s) it overlaps and correctly excludes the ones that don't overlap any gene.

**Wire to next stage:**
> Confirm Stage 9 reads this stage's ranked DUP/INV table for the report section specifically called out for Maria (duplications) and Alistair (inversions), and that the table is sorted/filterable by SV type so each of them can look at just their type of interest.

---

## Final section: integration, documentation, and demo

**Integration test:**
> Write an end-to-end integration test (e.g. `tests/test_integration.py` or a shell script `run_pipeline.sh`) that runs Stages 0 through 9 in sequence on synthetic data for a single chromosome, asserts each stage's expected output file exists and is non-empty, and asserts the final Stage 9 summary table has exactly one row per haploblock with all expected columns populated (no unexpected nulls). This should run in under a few minutes so it can be re-run after every change.

**Documentation:**
> Review `README.md`'s pipeline overview table against the actual `pipeline/` scripts as implemented, and update either the README or the code so stage names, inputs, and outputs match exactly. Add a `--help` docstring/argparse description to every stage script consistent with the README's one-line purpose for that stage.

**Demo / report:**
> Write `pipeline/stage9_integrate.py` that aggregates all prior stages' output tables into one per-haploblock summary (SV-type composition, boundary-enrichment flag, population-specificity density, cluster-correlation value, top DUP/INV hits) and renders it as a Jupyter notebook or a static HTML report (via Jinja2 + matplotlib/seaborn plots) — one figure per hypothesis (boundary enrichment histogram, per-block type-enrichment heatmap, population-correlation scatter, SV-based clustering scatter colored by both label sets, DUP/INV ranked table). If time allows, wrap the summary table in a minimal web view (e.g. a single-page app or a simple Flask/FastAPI endpoint) that lets a user query one haploblock ID and see its full SV profile, as a stretch extension of haploblocks.org itself.
