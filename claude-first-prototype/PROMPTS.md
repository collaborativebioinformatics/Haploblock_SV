# PROMPTS.md — Implementation handoff for the Haploblock_SV pipeline

Ready-to-use prompts for an agentic coding assistant, organized by pipeline stage and by day (see `README.md` for the full plan). Each prompt is self-contained enough to hand to a fresh agent session — it restates the relevant context rather than assuming prior conversation.

Shared context to paste into any session if it doesn't already have it: *"We're building a Python pipeline for a hackathon that studies structural variants (SVs: DEL/DUP/INV/INS) within 'haploblocks' — LD-defined haplotype-hash regions from data.haploblocks.org. SV calls come from dbVar study nstd152 (Chaisson et al. 2019), or from an arbitrary multi-sample SV VCF supplied to Stage 0 via `--vcf`. Population labels come from whatever is in `sample_metadata.tsv` (a `population` column) — never a hardcoded list of the five 1000G superpopulations. The pipeline is a sequence of independently-runnable, parameterized stages under `pipeline/`, each a script taking input paths and config values as CLI args — never hardcoded paths. Stages 0–2, 4, and 5 are implemented; Stage 3 was removed (numbering keeps 4–9); Stages 6–9 are specs. See README.md for the full stage list."*

---

## Day 1 (Aug 25) — Stages 0–2

### Stage 0: Data ingestion & harmonization

**Implement:**
> Write `pipeline/stage0_ingest.py`, a standalone script that: (1) obtains SV calls from one of, in priority order, `--vcf` (a standard multi-sample SV VCF with FORMAT/GT), `--sv-source` (a single prepared VCF), or an auto-fetch of the dbVar nstd152 call+region VCFs; (2) downloads or accepts a local path to haploblock boundaries from data.haploblocks.org; (3) obtains sample→population metadata from `--panel-source` if given, else derives it from the SV sample IDs; (4) confirms all coordinate inputs are on the same genome build (liftover to GRCh38 with `pyliftover` if a chain file is supplied, logging a warning otherwise); and (5) writes one shared `config.yaml` capturing genome build, thresholds (AF, boundary-distance N bp, size limits, `drop_imprecise`), seeds, and the three table paths — this config is the contract every later stage reads. All input paths/URLs are CLI args with sensible defaults, never hardcoded. Any real source that is slow, unreachable, or unparseable falls back to a small synthetic stand-in for just that input (logged clearly) so downstream stages are never blocked.
>
> For the `--vcf` path specifically: write `sv_calls.tsv` with **one column per sample holding the raw FORMAT/GT string** (`0|1`, `1/1`, `./.` — phasing and missingness preserved, unlike the 0/1/2 dosage the dbVar/synthetic paths write), and build `sample_metadata.tsv` from the VCF's sample list, with `population` = `UNKNOWN` unless `--panel-source` fills it in. `sample_metadata.tsv` always has `sample_id`, `population` (the fine-grained label Stages 4/6/7 group by), and `superpopulation` columns.
>
> `--drop-imprecise` sets `thresholds.drop_imprecise: true` in the emitted config; it is **off by default** because SV callers flag nearly every inversion IMPRECISE, so dropping imprecise calls removes all INV events.

**Test:**
> Write pytest tests for `stage0_ingest.py`. In synthetic-data mode (no network): `config.yaml` parses as YAML with all expected keys (genome build, thresholds incl. `drop_imprecise`, seeds, paths), the SV and haploblock tables are non-empty DataFrames with the expected columns, and `sample_metadata.tsv` has `sample_id` / `population` / `superpopulation`. In `--vcf` mode against a tiny checked-in example VCF: every per-sample cell in `sv_calls.tsv` is a GT string (contains `/` or `|`), the INV record is retained with `imprecise=True`, `config.yaml` has `drop_imprecise: false`, and `sample_metadata.tsv` is built from the VCF header with `population` all `UNKNOWN`.

**Wire to next stage:**
> Confirm `stage0_ingest.py`'s output file paths (SV table, haploblock table, sample metadata, config.yaml) match exactly what `pipeline/stage1_qc.py` expects to read — if Stage 1 doesn't exist yet, write down the expected input contract (file names, formats, columns) as a short comment block at the top of `stage0_ingest.py` so Stage 1 can be implemented against it independently.

### Stage 1: QC & normalization

**Implement:**
> Write `pipeline/stage1_qc.py`, a standalone script taking the raw SV table and haploblock table (paths from Stage 0's config.yaml) and: dropping SVs below/above configurable size thresholds (read from config.yaml; rows with unresolvable length are exempted, not dropped), coordinate sanity-checks (`start <= end`), deduplicating exact-duplicate calls (key on `chrom,start,end,sv_type,length` — INS share a ~1bp interval so `length` must be in the key), and validating the haploblock BED is sorted and non-overlapping per chromosome (raise a clear error listing offending rows, not silently continuing). Output a cleaned SV table, the haploblock table copied through unchanged, a QC report (counts before/after each step) as JSON, and this stage's own `config.yaml` with `paths` repointed so Stage 2 can chain via `--config`.
>
> **IMPRECISE calls are kept by default.** A `--drop-imprecise` flag (and `thresholds.drop_imprecise` in the config) enables the stricter filter, but it is off by default: callers flag nearly every inversion IMPRECISE, so dropping imprecise calls silently removes all INV events. True reference-based left-alignment is out of scope — accept `--reference-fasta` but make it an explicit logged no-op.

**Test:**
> Write pytest tests for `stage1_qc.py` chaining off Stage 0's synthetic output. Assert: the QC report's before/after counts are internally consistent; with defaults, imprecise calls are **retained**; with `--drop-imprecise`, all imprecise calls are removed and `qc_report.json`'s `dropped_imprecise` is > 0; size-filtered and duplicate rows are gone; the haploblock table is passed through byte-for-byte; and an overlapping haploblock table makes the stage exit non-zero with the offending rows named.

**Wire to next stage:**
> Update `pipeline/stage2_intersect.py` (or its docstring/contract if not yet written) to read Stage 1's cleaned SV table and cleaned haploblock table as its inputs, and confirm the coordinate convention (0-based BED-style, as fixed in Stage 1) is documented and consistent between the two scripts.

### Stage 2: SV × haploblock intersection

**Implement:**
> Write `pipeline/stage2_intersect.py` (pure Python/numpy — no bedtools/pybedtools) that intersects the QC'd SV table with the predefined haploblock BED. **For each SV, check it against every haploblock on its chromosome** with a plain vectorised interval comparison (an `[n_sv, n_block]` boolean grid) — no sorted-search window that could drop a match if block order were ever off. Classify by **overlap count** (half-open intervals, `sv.start < block.end and sv.end > block.start`): 0 blocks → `outside_block`, exactly 1 → `within_block`, ≥ 2 → `boundary_crossing` (the SV interval physically straddles the shared edge between contiguous blocks — its start is below one block's end while its end is above it / above the next block's start). **No proximity threshold**: an SV that lies wholly inside one block is `within_block` however close it sits to an edge. Add `haploblock_id` (comma-joined for `boundary_crossing`). Still re-validate the haploblock table as sorted/non-overlapping defensively, but do not let the classification depend on ordering. (`--boundary-distance-bp` / `thresholds.boundary_distance_bp` may still be accepted for backward compatibility but must be ignored.)
>
> Haploblocks are contiguous *within the span they cover* but do not reach the telomeres/centromere, so `outside_block` SVs are expected (telomeric calls). Log an `outside_block` breakdown — `before_first_block` / `after_last_block` / `in_inter_block_gap` / `no_blocks_on_chrom` — so a non-zero `in_inter_block_gap` (a real hole in the haploblock table) is visible rather than hidden. Output one annotated SV table (TSV) that is the shared input for downstream stages; take all paths as CLI args.

**Test:**
> Write pytest tests for `stage2_intersect.py` with a hand-crafted tiny table (two contiguous blocks + a lone block on another chromosome). Cover: an SV wholly inside one block → `within_block`; an SV a few hundred bp inside the *second* of two contiguous blocks, near the shared edge but not crossing it → `within_block` with that one block id (regression — a proximity rule used to mislabel this `boundary_crossing`); an SV straddling the shared edge → `boundary_crossing` with both ids; an SV overlapping nothing, and one on a chromosome with no blocks → `outside_block`. Also assert the `outside_block` breakdown log reports the right `before_first_block` / `after_last_block` / `in_inter_block_gap=0` / `no_blocks_on_chrom` counts on a table whose blocks start well after coordinate 0.

**Wire to next stage:**
> Confirm the annotated SV table's schema (SV id, type, chrom, start, end, haploblock_id, position_class, per-sample genotype columns) is exactly what Stage 4 (AF classification) needs — document it as a comment block at the top of `stage2_intersect.py`, including that per-sample cells may be GT strings (from `--vcf`) or dosage ints.

---

## Day 2 (Aug 26) — Stages 4–5

> Stage 3 (boundary enrichment test) was removed from the pipeline — see README.md "Descoped / future steps". Stage numbers below are unchanged (4–9) so these prompts and the README stay aligned.

### Stage 4: Common vs. population-specific SV classification

**Implement:**
> Write `pipeline/stage4_classify_af.py` that reads Stage 2's annotated SV table (per-sample genotype columns — handle both GT strings like `0|1` and 0/1/2 dosage ints) and `sample_metadata.tsv`, and computes allele frequency per SV **within each population present in `sample_metadata.tsv`'s `population` column** — do not hardcode the five 1000G superpopulations, and do not use the data.haploblocks.org clusters here (that would make Stage 6 circular). For each `(sv_id, sv_type, haploblock_id)` emit the per-population AF and a category: `common` (AF ≥ a configurable threshold, default 0.05, in ≥ 2 populations), `specific_to_population` (AF ≥ threshold in exactly one population, near-zero elsewhere — record which population), or `other` (neither — typically too few samples in the relevant populations to tell; log the fraction here). Write the result as a tidy table.

**Test:**
> Write a pytest test with a synthetic 3-SV, 2-population genotype matrix: one SV common to both populations, one private to population A, one rare/absent everywhere. Assert the category (`common`, `specific_to_population` with the correct population, `other`) matches for all three at the default threshold, and that the per-population AF columns are correct.

**Wire to next stage:**
> Stage 6 (population-cluster correlation) and Stage 8 (DUP/INV overlay) both consume `sv_af_classification.tsv` — the `sv_category` column (value `specific_to_population`), the `specific_to_population` column (which population), `sv_type`, and `haploblock_id`. Both must first collapse the tidy table to one row per `sv_id` (these fields are constant across a SV's per-population rows). This table carries **no** SNV-cluster / hash column and **no** per-sample genotypes: Stage 6 reads its cluster labels only from its own `--clusters` input (keeps the correlation non-circular), and Stage 8 joins back to Stage 2's `sv_calls.tsv` on `sv_id` for the genotypes it needs to count carriers/recurrence.

**Implemented:** `pipeline/stage4_classify_af.py`. Output `sv_af_classification.tsv` is tidy/long — one row per (SV × population) with `af`, `n_called`, `pop_has_data`, and per-SV `sv_category` / `specific_to_population` / `other_reason` (`absent_or_rare` / `one_pop_high_plus_intermediate_elsewhere` / `insufficient_population_data`). Flags: `--af-threshold` (default = config's `af_common_threshold`), `--absent-af-threshold` (0.01), `--min-samples-per-pop` (2). Example: `example_data/stage4_example/`. `haploblock_id` is passed through from Stage 2 verbatim (comma-joined ids for `boundary_crossing` SVs are not split).

### Stage 5: Per-haploblock SV-type enrichment

**Implement:**
> Write `pipeline/stage5_type_enrichment.py` that builds a per-haploblock × SV-type count matrix from Stage 2's output. For each SV type, calculate the overall SV rate across all haploblocks as the total number of SVs of that type divided by the total haploblock length. For each haploblock, calculate the expected number of SVs of that type based on its length and the overall SV rate. Use a Poisson test to determine whether the observed SV count significantly deviates from the expected count. Apply Benjamini-Hochberg FDR correction across all haploblock × SV-type tests. Flag combinations with an FDR-adjusted q-value below 0.05. Output a results table: `haploblock_id, sv_type, observed_count, expected_count, p_value, q_value, flagged`.

**Test:**
> Write a pytest test with synthetic data: several haploblocks of varying length with SV counts proportional to length (should NOT be flagged after the length adjustment), plus one haploblock with an artificially inflated count for one SV type (should be flagged after FDR correction). Assert the inflated block×type is flagged and the proportional ones are not, at a fixed seed.

**Wire to next stage:**
> Stage 9 (integration) reads `sv_type_enrichment.tsv` (path via Stage 5's `config.yaml` → `paths.sv_type_enrichment`) for the "haploblocks prone / resistant to a particular SV type" section. The table has no direction column — Stage 9 derives it from the sign of `observed_count − expected_count` on rows where `flagged` is true: `observed > expected` → *prone*, `observed < expected` → *resistant*. The per-block type-enrichment heatmap uses `log2((observed_count + 0.5) / (expected_count + 0.5))` per (haploblock, sv_type) cell, with flagged cells marked.

**Implemented:** `pipeline/stage5_type_enrichment.py`. `outside_block` SVs dropped; `boundary_crossing` SVs with a comma-joined `haploblock_id` counted once per spanned block. Exact two-sided Poisson test (`2·min(P(X≤obs), P(X≥obs))`, capped at 1); BH-FDR via `scipy.stats.false_discovery_control` (statsmodels fallback) across the **full** haploblocks × observed-types grid, zeros included. `--q-threshold` (default = config's `thresholds.q_threshold`, else 0.05). Output `sv_type_enrichment.tsv` = `haploblock_id, sv_type, observed_count, expected_count, p_value, q_value, flagged`, sorted by `q_value`. Example + generator: `example_data/stage5_example/`.

---

## Day 3 (Aug 27) — Stages 6–8 (parallelizable)

### Stage 6: Population-cluster correlation

**Implement:**
> Write `pipeline/stage6_cluster_correlation.py` that evaluates how well Stage 4's population-specific SV calls correspond to the SNV-based population clusters already published for these haploblocks. Read Stage 4's `sv_af_classification.tsv` and **collapse it to one row per `sv_id`** (`sv_category`, `specific_to_population`, `sv_type`, `haploblock_id` are constant across a SV's per-population rows). Per haploblock, compute a population-specific-SV density = count of SVs with `sv_category == "specific_to_population"` / total SVs (or / block length); for `boundary_crossing` SVs whose `haploblock_id` is comma-joined, either explode to each block or restrict this stage to `position_class == "within_block"` — state which.
>
> The SNV-based cluster labels come **only** from a `--clusters` flag that builds a per-haploblock cluster table from a data.haploblocks.org clusters file (local path or URL; ideally fetched once in Stage 0 and pointed at here). Never read a cluster/hash column from Stage 4's output — it has none, by design, and that is what keeps this correlation non-circular. Compare the two per-block quantities with a Spearman correlation and a permutation p-value (shuffle block labels, default 1000 permutations, seed from config.yaml). Output the correlation coefficient, p-value, and a per-block table for plotting.

**Test:**
> Write a pytest test with synthetic data where population-specific-SV density is a noisy monotonic function of a synthetic per-block cluster score; assert the Spearman correlation is positive and significant (p < 0.05) at a fixed seed, and that a shuffled-label negative control on the same data is not significant.

**Wire to next stage:**
> Confirm Stage 9 reads this stage's correlation coefficient/p-value and per-block table to render the population-correlation section of the report.

### Stage 7: SV-based population structure reconstruction

**Implement:**
> Write `pipeline/stage7_sv_clustering.py` that builds a per-sample × per-haploblock SV presence/absence (or dosage) matrix from Stage 2's output and genotype columns, runs PCA (scikit-learn) and UMAP (umap-learn, fixed `random_state` from config.yaml), and clusters the reduced embedding (e.g. k-means or HDBSCAN). **Output only: PNG PCA and UMAP plots (points colored by the `population` labels from `sample_metadata.tsv`), plus a cluster-assignments table (`sample_id, pca_cluster, umap_cluster`).** Formal cluster-agreement statistics (ARI vs. populations and vs. the haploblocks.org hash clusters) are a follow-up — see README.md "Descoped / future steps".

**Test:**
> Write a pytest test with a small synthetic SV matrix where two groups of samples have clearly distinct SV profiles; assert the clustering recovers the two groups (each synthetic group maps predominantly to one cluster label) at a fixed seed, and that the expected PNG files and the cluster-assignments table are written.

**Wire to next stage:**
> Confirm Stage 9 can read the cluster-assignments table and embed the PNGs in the final report.

### Stage 8: Duplication/inversion gene overlay

**Implement:**
> Write `pipeline/stage8_dup_inv_overlay.py` that starts from Stage 4's `sv_af_classification.tsv` collapsed to one row per `sv_id`, keeps `sv_type in {DUP, INV}`, and reads the `sv_category` / `specific_to_population` columns from there. Recurrence (present in ≥ 2 samples) is **not** in Stage 4's output — join back to Stage 2's `sv_calls.tsv` on `sv_id` for the per-sample genotype columns and count carriers there. Flag recurrent calls and `sv_category == "specific_to_population"` calls, and overlap them against a gene annotation GTF (via `pyranges`/`gtfparse`) to report which genes each overlaps. Output a ranked table (by carrier count, then by population specificity) of "stand-out" DUP/INV calls with overlapping gene names. (A selection-scan overlay is a follow-up — see README.md "Descoped / future steps".)

**Test:**
> Write a pytest test with a synthetic DUP/INV table and a tiny synthetic GTF (2-3 genes at known coordinates); assert the gene-overlap output correctly attributes each synthetic SV to the gene(s) it overlaps and excludes the ones that don't overlap any gene.

**Wire to next stage:**
> Confirm Stage 9 reads this stage's ranked DUP/INV table for the report section called out for Maria (duplications) and Alistair (inversions), sorted/filterable by SV type.

---

## Final section: integration, documentation, and demo

**Integration test:**
> Write an end-to-end integration test (e.g. `tests/test_integration.py` or a shell script `run_pipeline.sh`) that runs the retained stages (0, 1, 2, then 4–9) in sequence on synthetic data for a single chromosome, asserts each stage's expected output file exists and is non-empty, and asserts the final Stage 9 summary table has exactly one row per haploblock with all expected columns populated (no unexpected nulls). This should run in under a few minutes so it can be re-run after every change.

**Documentation:**
> Review `README.md`'s pipeline overview table against the actual `pipeline/` scripts as implemented, and update either the README or the code so stage names, inputs, and outputs match exactly. Add a `--help` docstring/argparse description to every stage script consistent with the README's one-line purpose for that stage.

**Demo / report:**
> Write `pipeline/stage9_integrate.py` that aggregates the prior stages' output tables into one per-haploblock summary and renders it as a Jupyter notebook or a static HTML report (Jinja2 + matplotlib/seaborn). Inputs, each located via its stage's `config.yaml` `paths.*` entry:
> - Stage 5 `sv_type_enrichment.tsv` → the "haploblocks prone / resistant to a particular SV type" section: filter `flagged == True`, split by sign of `observed_count − expected_count` (`>` = prone, `<` = resistant), list `haploblock_id, sv_type, observed_count, expected_count, q_value`; and a per-block × SV-type heatmap of `log2((observed_count + 0.5)/(expected_count + 0.5))` with flagged cells marked.
> - Stage 4 `sv_af_classification.tsv` → per-haploblock population-specificity density (share of SVs with `sv_category == "specific_to_population"`).
> - Stage 6 correlation coefficient / p-value + per-block table → population-correlation scatter.
> - Stage 7 cluster-assignments table + PNGs → SV-based PCA/UMAP figure.
> - Stage 8 ranked DUP/INV table → the DUP/INV hits section.
> Emit one per-haploblock summary row per block with all of the above columns populated.
