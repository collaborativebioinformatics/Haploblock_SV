# Haploblock_SV — Structural Variants Within Haploblocks

Prototype pipeline plan for the Structural Variants Hackathon at Baylor College of Medicine, August 25–28, 2026. Stages 0–2 are implemented and run end-to-end on the example data (see the "Prototype" sections below and `PROMPTS.md` for the per-stage handoff prompts). Stages 4–9 are still prompt specs pending biological review. Stage 3 (boundary enrichment) has been dropped from the pipeline — see "Descoped / future steps".

## Introduction

Haploblocks — genomic regions of conserved haplotype structure identified by the [haploblocks.org](https://haploblocks.org) / [data.haploblocks.org](https://data.haploblocks.org) projects via genomic hashing — capture a layer of population structure that is complementary to single-SNP analyses, but how structural variants (SVs) are organized relative to these blocks is unexplored. Our core research question is whether SV type (deletion, duplication, inversion, insertion), location (within-block vs. boundary-crossing), and population specificity are non-randomly distributed across haploblocks, and whether that distribution tracks the population clusters haploblocks.org already derives from haplotype hashes.

Primary data is the 1000 Genomes long-read/haplotype-resolved SV callset at dbVar (Study ID **nstd152**, Chaisson et al. 2019); Stage 0 can also ingest an arbitrary multi-sample SV VCF via `--vcf` (keeping the raw genotypes), joined against haploblock coordinates and metadata from data.haploblocks.org. Population labels come from `sample_metadata.tsv` — whatever populations that file actually contains, not a hardcoded set of five superpopulations. The pipeline intersects SVs with haploblocks (Stage 2), classifies each SV by population specificity (Stage 4), then runs per-haploblock SV-type enrichment (Stage 5), a correlation against the pre-existing SNV-based clusters (Stage 6), an SV-based PCA/UMAP (Stage 7), and a DUP/INV gene overlay (Stage 8). What's novel is treating haploblocks as the unit of SV analysis rather than genes or fixed windows.

The expected deliverable is a reusable, parameterized pipeline plus a per-haploblock summary report (SV-type composition, population specificity, cluster correlation, and a ranked duplication/inversion list). Known risks going in: nstd152 covers only 9 samples across 3 superpopulations (population-level analyses are underpowered without a wider cohort via `--vcf`), per-block SV counts may be underpowered for rare SV types in small blocks, and validating findings against an independent long-read callset was scoped out as too data-wrangling-heavy for 3 days. Heavier statistical and visualization ideas (permutation-based boundary enrichment, negative-binomial offsets, cluster-agreement indices, selection-scan overlays, a web view) are parked in "Descoped / future steps".

## Pipeline overview

| # | Stage | Purpose | Status |
|---|---|---|---|
| 0 | Data ingestion & harmonization | Fetch/parse the nstd152 VCF (or ingest any multi-sample SV VCF via `--vcf`, keeping the raw GT), haploblock boundaries from data.haploblocks.org, and sample→population metadata; write one shared config (thresholds, genome build, seeds) | Implemented |
| 1 | QC & normalization | Filter SVs on size and coordinate sanity, dedup, validate the haploblock BED is sorted/non-overlapping. IMPRECISE calls are **kept by default** (dropping them removes every inversion); `--drop-imprecise` opts back in | Implemented |
| 2 | SV × haploblock intersection | Check every SV against every haploblock on its chromosome; classify by overlap count — 1 block = `within_block`, ≥ 2 = `boundary_crossing` (spans a shared edge), 0 = `outside_block`. No proximity threshold. Blocks are contiguous within their span but do not reach the telomeres, so `outside_block` = telomeric/centromeric SVs (logged with a before/after/in-gap breakdown) | Implemented |
| ~~3~~ | ~~Boundary enrichment test~~ | Removed from the pipeline (see "Descoped / future steps"). Stage numbers 4–9 are kept as-is for continuity with the prompts | Removed |
| 4 | Common vs. population-specific SV classification | Allele frequency per population from `sample_metadata.tsv` (whatever populations it holds — no fixed AFR/AMR/EAS/EUR/SAS list). Per `sv_id × sv_type × haploblock_id`: the per-population AF and a category of `common`, `specific_to_population`, or `other` (too little data) | Proposed |
| 5 | Per-haploblock SV-type enrichment | Per-haploblock × SV-type count matrix; overall per-type rate = total SVs of that type / total haploblock length; expected count per block from its length; Poisson test observed vs. expected; Benjamini-Hochberg FDR across all block × type tests; flag q < 0.05 | Proposed |
| 6 | Population-cluster correlation | Compare Stage 4's population-specific SV patterns against the predefined SNV-based haploblock clusters; `--clusters` builds the cluster table from a data.haploblocks.org clusters file | Proposed |
| 7 | SV-based population structure reconstruction | Per-sample × per-haploblock SV matrix → PCA/UMAP; output **PNG plots and cluster assignments only**, colored by the populations in `sample_metadata.tsv` | Proposed |
| 8 | Duplication/inversion gene overlay | Gene overlap for recurrent / population-specific DUP and INV calls | Proposed |
| 9 | Integration & report | Aggregate stage outputs into a per-haploblock summary table plus plots | Proposed |

Stages 6, 7, and 8 depend only on Stages 1/4 respectively and can run in parallel once QC is done — good for splitting across the team.

## Development plan

**Day 1 (Aug 25).** Stand up Stages 0–2. Timebox Stage 0 hard (half a day) with a synthetic-data fallback if the dbVar VCF or the data.haploblocks.org export is slow to parse, so nothing downstream is blocked. Get the SV × haploblock intersection (Stage 2) producing an annotated table on real or synthetic data by end of day.

**Day 2 (Aug 26).** Stages 4–5: population-specific classification and the length-adjusted per-block SV-type enrichment scan (Poisson + BH-FDR). This is the statistical core — get FDR-corrected results on the real dataset, even if plots/report are rough.

**Day 3 (Aug 27).** Stages 6–8 in parallel across the team (population-cluster correlation, SV-based PCA/UMAP, gene overlay for Maria/Alistair's duplication/inversion question). Start Stage 9 integration in the afternoon so there's a working summary table before the final day.

**Day 4 / buffer (Aug 28).** Finish Stage 9 (report + plots), run the end-to-end integration test below, prepare the demo/writeup. Treat this as buffer, not new-development time.

## Testing steps

- **Per-stage sanity checks on a small slice:** run the stages on one chromosome (e.g., chr21) or the synthetic haploblock+SV set before running genome-wide; confirm intersection counts (Stage 2), classification counts (Stage 4), and enrichment q-value distributions (Stage 5) look sane (e.g., q-values roughly uniform under a shuffled-label negative control).
- **Negative control for Stage 5:** re-run the enrichment test on label-shuffled data — it should show no flagged block × type combinations; a positive result on shuffled data indicates a pipeline bug, not a biological signal.
- **End-to-end integration run:** execute the retained stages on the full nstd152 + haploblocks.org dataset (or the `--vcf` example) for one chromosome, confirm the pipeline completes without manual intervention and the final summary table has one row per haploblock with all expected columns populated.
- **Validation against known structure:** in Stage 7, confirm the SV-based PCA/UMAP separates the populations present in `sample_metadata.tsv` before reading anything into finer structure — if it can't, the SV genotype matrix itself is broken.
- **Reproducibility check:** re-run the pipeline twice with the same config (fixed seeds) and confirm identical output.

## Non-pip dependencies

- `bcftools` / `samtools` (CLI) — optional; only needed if reference-based left-alignment is added to Stage 1 later (`bcftools norm -f`). The implemented Stage 1 does not call them.
- R is not required; the plan above is Python-only. If a team member prefers R for Stage 5/6 stats, that's a drop-in alternative to `scipy`/`statsmodels`, not a pipeline dependency.

## Prototype: Stage 0 (data ingestion)

`pipeline/stage0_ingest.py` is implemented and working — it's the current state of the prototype. **By default it fetches real data, genome-wide**: dbVar's nstd152 call+region VCFs, and haploblock boundaries for all 23 chromosomes (1-22 + X) from data.haploblocks.org. Any source that fails to download or parse falls back to synthetic data for just that input, so later stages are never blocked.

**Setup.** The `pyEnv_SVhack2026` virtualenv (`~/pyenvs/pyEnv_SVhack2026`) already has everything in `requirements.txt` installed. Elsewhere: `pip install -r requirements.txt` (plus `bcftools`/`samtools` for later stages).

**Run command** (no arguments needed — fetches real data genome-wide by default, ~14s):
```
~/pyenvs/pyEnv_SVhack2026/bin/python pipeline/stage0_ingest.py
```
This writes `example_data/sv_calls.tsv`, `example_data/haploblocks.tsv`, `example_data/sample_metadata.tsv`, and `example_data/config.yaml`, plus the four raw dbVar files (`nstd152.GRCh38.variant_{call,region}.vcf.gz[.tbi]`). Haploblock boundaries come from one small per-chromosome file each — `data.haploblocks.org` publishes `<chrom>/<chrom>_haploblock_boundaries_<chrom>.tsv` (one row per block) and usually `<chrom>/<chrom>_haploblock_hashes.tsv` (one representative hash per block) — rather than the thousands of individual per-block `chrN_cluster_hashes_<start>-<end>.tsv` files also present in the same directories, which don't scale to a genome-wide default (one HTTP request per block; an earlier version of this script that used those took 2-3 minutes for a 50-blocks-per-chromosome cap alone). Useful flags:
- `--vcf PATH` ingests a standard multi-sample SV VCF (FORMAT/GT) as the SV source. It is the highest-priority SV source (over `--sv-source` and the dbVar auto-fetch): `sv_calls.tsv` gets one column per sample holding the **raw GT string** (`0|1`, `1/1`, `./.` — phasing and missingness preserved, unlike the 0/1/2 dosage the other paths write), and `sample_metadata.tsv` is built from the VCF's sample list with `population` = `UNKNOWN` unless `--panel-source` is also given. Example: `--vcf example_data/example_cohort.vcf --haploblock-source example_data/haploblocks.tsv`.
- `--haploblock-chroms chr21,chr22` restricts to specific chromosomes instead of the `all` default.
- `--skip-dbvar-download` / `--skip-haploblock-download` force synthetic data for that input, for fast offline iteration.
- `--sv-source` / `--haploblock-source` / `--panel-source` point at a single already-prepared local file or URL instead (bypasses the auto-fetchers above).
- `--drop-imprecise` sets `thresholds.drop_imprecise: true` in the config so Stage 1 drops IMPRECISE calls. **Off by default** — callers flag nearly every inversion IMPRECISE, so dropping them silently removes all INV events.

**Expected output** — log lines naming which of the three inputs are real vs. synthetic and key data-quality warnings, e.g.:
```
INFO: chr1: fetched 3049 haploblock boundaries
...
INFO: chrX: fetched 846 haploblock boundaries
INFO: Fetched 39113 haploblock boundaries across 23 chromosome(s)
INFO: dbVar call file: 102153 usable (sample, call) observation(s) across 9 sample(s)
INFO: Using real SV calls (50601 records) from ... (nstd152)
WARNING: No --panel-source given: replaced the placeholder synthetic sample panel with 9 sample(s)
  inferred from the real dbVar SV calls (hardcoded HGSVC trio lookup). This covers only
  ['AFR', 'AMR', 'EAS'] -- population-specific classification cannot be evaluated for any
  superpopulation not listed.
INFO: Wrote 50601 SVs, 39113 haploblocks, 9 samples to .../example_data
```
`example_data/config.yaml` is the contract every later stage reads (genome build, AF/boundary thresholds, seeds, and the three table paths) — inspect it to confirm it parses and has all five top-level keys (`genome_build`, `data_sources`, `thresholds`, `seeds`, `paths`).

**Important data-availability caveat, not a bug:** nstd152's `variant_call.vcf.gz` currently covers exactly **9 samples from 3 HGSVC trios** (CHS: HG00512/13/14, PUR: HG00731/32/33, YRI: NA19238/39/40) — three populations only. Without a real `--panel-source`, Stage 0 infers labels for just these 9 via hardcoded lookups in `stage0_ingest.py` (`sample_metadata.tsv` gets `population` = CHS/PUR/YRI and `superpopulation` = EAS/AMR/AFR). Stage 4/6/7's population-level analyses are severely underpowered on this alone (n=3 per population) — supply a wider cohort with `--vcf` (plus `--panel-source` for its population labels) to make them meaningful. nstd152 itself only has calls for these 9 individuals, so a panel alone does not add samples.

**Known limitations, deliberately out of scope for this prototype:** no retry/resume on failed downloads; `n_clusters` (the number of distinct haplotype clusters in a block) is always `NaN` for real haploblock data since it needs the per-block cluster-enumeration files this script deliberately avoids fetching at genome scale — only `hash_length` (derived from the one representative hash per block) is recovered; `cluster_diff_score` (an actual population-differentiation metric) is left as `NaN` for real haploblock data — computing it needs population labels, which belongs in a later stage (Stage 6), not ingestion; genome build is taken from explicit `--*-genome-build`/`--dbvar-build` flags rather than auto-detected from file headers.

`sv_calls.tsv` also carries `imprecise` (bool) and `length` (float, nullable) columns — dbVar never sets VCF QUAL/FILTER (always `"."` in nstd152), so `imprecise` is the closest thing this source has to a confidence flag, and INS records need `length` because their `end - start` is always ~1bp (an insertion's reference position, not its size). Per-sample columns hold a 0/1/2 dosage int for the dbVar and synthetic paths, but the **raw GT string** (`0|1`, `1/1`, `./.`) for the `--vcf` path. See the "Stage 1" comment block at the top of `stage0_ingest.py` for the full column contract.

## Prototype: Stage 1 (QC & normalization)

`pipeline/stage1_qc.py` is implemented and working. It reads Stage 0's `config.yaml`, applies a size filter (`min_sv_length`/`max_sv_length` from config, INS rows with unresolvable length are exempted rather than dropped), coordinate sanity-checks (`start <= end`), and dedup, then validates the haploblock table is sorted and non-overlapping per chromosome — raising with the specific offending rows listed, not silently continuing, if it isn't.

**IMPRECISE calls are kept by default.** Most SV callers flag essentially every inversion (and many insertions) IMPRECISE because their breakpoints are fuzzy — so dropping imprecise calls silently removes *all* INV events, which is wrong for a study of SV-type distribution. Pass `--drop-imprecise` (or set `thresholds.drop_imprecise: true` in the config) to opt into the stricter filter.

**Dedup note, learned from the real data:** deduping on `(chrom, start, end, sv_type)` alone is wrong — 2439/50601 real rows collide that way, almost all INS, because dbVar gives every insertion the same ~1bp reference interval regardless of its actual size, so two textually "identical" rows are usually two distinct real insertions at the same locus. Adding `length` to the key resolves nearly all of that; the ~600 groups that still collide even with `length` included turn out to be the same event reported twice at different confidence (one precise, one imprecise) — ties are broken by keeping the precise record.

**Run command** (chains off Stage 0's output):
```
~/pyenvs/pyEnv_SVhack2026/bin/python pipeline/stage1_qc.py --config example_data/config.yaml --out-dir stage1_output
```
Writes `stage1_output/sv_calls.tsv`, `haploblocks.tsv`, `sample_metadata.tsv` (copied through), `qc_report.json` (before/after counts per filter step), and its own `config.yaml` (same schema as Stage 0's, `paths` repointed at these cleaned files) so Stage 2 can chain via `--config` the same way. `--min-sv-length`/`--max-sv-length`/`--drop-imprecise` override Stage 0's config values without re-running ingestion.

**Expected output** on the real genome-wide nstd152 data pulled by Stage 0's default run (imprecise kept):
```
INFO: Confidence filter: dropped 0 imprecise row(s) (drop_imprecise=False)
INFO: Size filter [50, 5000000]bp: dropped 31 row(s), exempted 2 with unknown length
INFO: Dedup: removed 612 exact-duplicate row(s)
INFO: Haploblock validation passed: 39113 block(s), sorted and non-overlapping per chromosome
INFO: Wrote 49958/50601 SVs and 39113 haploblocks to .../stage1_output
```

**Known limitations, deliberately out of scope:** true left-alignment (VCF-style indel normalization) needs the reference genome; `--reference-fasta` is accepted but is an explicit no-op (logged loudly) — `bcftools norm -f` is the standard tool for this if/when a reference is added.

## Prototype: Stage 2 (SV × haploblock intersection)

`pipeline/stage2_intersect.py` is implemented and working — **pure Python/numpy, no bedtools/pybedtools**. For each chromosome it compares *every* SV against *every* haploblock on that chromosome with a vectorised interval comparison (an `[n_sv, n_block]` boolean grid). There is no sorted-search window that could drop a match if block order were ever off; blocks are still re-validated as sorted and non-overlapping, but the classification no longer depends on that. Genome-wide, ~50k SVs against 39,113 haploblocks classify in ~2 s.

Each SV gets `position_class` and `haploblock_id` from a pure **overlap count**: overlaps 0 blocks → `outside_block`; exactly 1 → `within_block`; ≥ 2 → `boundary_crossing` (the SV interval physically straddles the shared edge(s) between contiguous blocks, so `haploblock_id` lists all of them). There is **no proximity threshold** — an SV sitting 600 bp inside one block, close to an edge but not crossing it, is `within_block`, not `boundary_crossing`. (`--boundary-distance-bp` / `thresholds.boundary_distance_bp` are still accepted so older configs don't error, but are ignored.)

**On `outside_block` — not a bug:** real data.haploblocks.org blocks are contiguous *within the span they cover* (no inter-block gaps, confirmed genome-wide), but that span does not reach the telomeres/centromere (e.g. chr21's blocks span ~14.2–46.2 Mb). SVs before the first block or after the last block on their chromosome are legitimately `outside_block`. Stage 2 logs a breakdown — `before_first_block` / `after_last_block` / `in_inter_block_gap` / `no_blocks_on_chrom` — so a genuine gap in the haploblock table (a non-zero `in_inter_block_gap`) would stand out instead of hiding among the expected telomeric calls.

**Run command** (chains off Stage 1's output):
```
~/pyenvs/pyEnv_SVhack2026/bin/python pipeline/stage2_intersect.py --config stage1_output/config.yaml --out-dir stage2_output
```

**Expected output** on the real genome-wide haploblocks + nstd152 data (imprecise kept in Stage 1):
```
INFO: Note: boundary_distance_bp is no longer used -- an SV is boundary_crossing iff it overlaps >=2 haploblocks
INFO: Position classification: within_block=45393 boundary_crossing=1519 outside_block=3046
INFO: 1519 SV(s) overlap more than one haploblock (span a shared edge)
INFO: outside_block breakdown: before_first_block=1692 after_last_block=1339 in_inter_block_gap=0 no_blocks_on_chrom=15
```
`in_inter_block_gap=0` is the check that matters: every `outside_block` SV is a telomeric call before the first block or after the last, or (15 of them) on a chromosome with no haploblocks at all (chrY) — none fell in a gap *between* blocks, confirming the haploblock table is contiguous and the matching has no ordering bug. Only 1,519 SVs genuinely straddle a block edge; the rest sit cleanly inside one block.

## Tests

`tests/test_stage0_ingest.py`, `tests/test_stage1_qc.py`, and `tests/test_stage2_intersect.py` run all three stages end-to-end on synthetic (offline) data — plus a `--vcf` path check, the imprecise-kept-by-default / `--drop-imprecise` behaviour, hand-crafted small-table cases for Stage 2's exact classification rules, and the `outside_block` breakdown log. Run with:
```
~/pyenvs/pyEnv_SVhack2026/bin/python -m pytest tests/ -v
```

## Descoped / future steps

Moved out of the hackathon pipeline to keep the prototype focused; each is a self-contained follow-up study:

- **Stage 3 — boundary enrichment test.** Permutation null (circular rotation of SV positions per chromosome) for the count of boundary-crossing SVs, plus a nearest-neighbor / Ripley's-K spacing check per block. Removed because the descriptive `boundary_crossing` counts from Stage 2 answer the core question for now, and a correct permutation null needs a per-chromosome covered-span model.
- **Richer Stage 5 model.** Negative-binomial (not just Poisson) with a dispersion test, `log(block_length)` offset *plus* SNP-density covariate, and a configurable minimum-count exclusion with an "underpowered fraction" report. The implemented Stage 5 spec keeps just the length-adjusted Poisson + BH-FDR.
- **Haploblock complexity vs. SV type.** Test whether the number of haplotype clusters in a block relates to the SV types it carries — i.e. whether DEL/INS/INV/DUP enrichment tracks block complexity. Needs the per-block cluster-enumeration files Stage 0 deliberately does not fetch at genome scale.
- **Cluster-agreement statistics for Stage 7.** Adjusted Rand Index (and similar) between the SV-based clusters and both the `sample_metadata.tsv` populations and the haploblocks.org hash clusters. Stage 7's core output is just the PCA/UMAP PNGs and cluster assignments.
- **Selection-scan overlay (Stage 8b).** Cross-reference recurrent DUP/INV calls against a selection-scan region file (e.g. PopHumanScan BED), adding a `near_selection_signal` column.
- **Web view.** A single-page app / small FastAPI endpoint to query one haploblock ID and see its full SV profile, as an extension of haploblocks.org.
