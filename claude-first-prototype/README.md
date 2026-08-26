# Haploblock_SV — Structural Variants Within Haploblocks

Prototype pipeline plan for the Structural Variants Hackathon at Baylor College of Medicine, August 25–28, 2026. This document is the output of a brainstorming/planning pass (see `PROMPTS.md` for implementation handoff prompts); no pipeline code has been written yet.

## Introduction

Haploblocks — genomic regions of conserved haplotype structure identified by the [haploblocks.org](https://haploblocks.org) / [data.haploblocks.org](https://data.haploblocks.org) projects via genomic hashing — capture a layer of population structure that is complementary to single-SNP analyses, but how structural variants (SVs) are organized relative to these blocks is unexplored. Our core research question is whether SV type (deletion, duplication, inversion, insertion), location (within-block vs. boundary-crossing), and population specificity are non-randomly distributed across haploblocks, and whether that distribution tracks the population clusters haploblocks.org already derives from haplotype hashes. The interdisciplinary angle draws spatial-statistics methods (permutation nulls, nearest-neighbor spacing tests borrowed from recombination-hotspot analysis) into the boundary-enrichment question, and dimensionality-reduction/cluster-agreement metrics (PCA/UMAP + Adjusted Rand Index) into the population-structure question, so that "correlates with population clusters" becomes a number rather than a visual impression. Primary data is the 1000 Genomes long-read/haplotype-resolved SV callset at dbVar (Study ID **nstd152**, Chaisson et al. 2019), joined against haploblock coordinates, metadata, and population-cluster labels from data.haploblocks.org and standard 1000 Genomes superpopulation assignments. The approach intersects SVs with haploblocks, classifies each SV by position (within-block, boundary-crossing) and population specificity (common vs. population-specific), then runs enrichment, regression, and clustering analyses on top of that shared annotated table. What's novel is treating haploblocks as the unit of SV analysis rather than genes or fixed windows, and testing — rather than assuming — whether SV-derived population structure agrees with the project's own existing hash-based clusters. The expected deliverable is a reusable, parameterized pipeline plus a per-haploblock summary report (SV-type composition, boundary flags, population specificity, cluster correlation, and a ranked duplication/inversion list of interest to the team). Known risks going in: the exact schema of the data.haploblocks.org cluster export is unverified until Stage 0 runs, per-block SV counts may be underpowered for rare SV types in small blocks, and validating findings against an independent long-read (ONT/Sniffles2) callset was scoped out of the hackathon as too data-wrangling-heavy for 3 days.

## Pipeline overview

| # | Stage | Purpose | Addresses |
|---|---|---|---|
| 0 | Data ingestion & harmonization | Fetch/parse nstd152 VCF, haploblock BED + metadata + cluster labels, 1KGP population metadata; write one shared config (thresholds, genome build, random seeds) | Foundation |
| 1 | QC & normalization | Filter SVs (PASS, size/confidence thresholds), normalize breakpoints, sanity-check haploblock BED | Foundation |
| 2 | SV × haploblock intersection | Classify each SV as within-block, boundary-crossing (breakpoint within N bp of an edge), or outside | H1 |
| 3 | Boundary enrichment test | Permutation test (circular rotation of SV positions per chromosome) + nearest-neighbor/Ripley's-K spacing check | H1, H7 |
| 4 | Common vs. population-specific SV classification | Per-superpopulation allele frequency → common (≥2 superpops) vs. population-specific (private), using standard AFR/AMR/EAS/EUR/SAS labels, kept independent of haploblocks.org clusters | H3, H5 |
| 5 | Per-haploblock SV-type enrichment | Poisson/negative-binomial regression, block length + SNP density as offset, per-block-per-type deviation, BH-FDR corrected, minimum-count threshold | H2, H4 |
| 6 | Population-cluster correlation | Correlate per-block population-specific-SV density against data.haploblocks.org's cluster differentiation metric | H3 |
| 7 | SV-based population structure reconstruction | Per-sample/per-haploblock SV matrix → PCA/UMAP → Adjusted Rand Index vs. superpopulation labels and vs. haploblocks.org clusters | H6 |
| 8 | Duplication/inversion gene & selection overlay | Gene overlap (core) for recurrent/population-specific DUP/INV; selection-scan cross-reference (8b, optional/stretch) | H5 |
| 9 | Integration & report | Aggregate all stage outputs into a per-haploblock summary table, plots, and (stretch) a simple web view | Deliverable |

Stages 6, 7, and 8 depend only on Stages 1/4 respectively and can run in parallel once QC is done — good for splitting across the team.

## Development plan

**Day 1 (Aug 25).** Stand up Stages 0–2. Timebox Stage 0 hard (half a day) with a synthetic-data fallback if the dbVar VCF or the data.haploblocks.org export is slow to parse, so nothing downstream is blocked. Get the SV × haploblock intersection (Stage 2) producing an annotated table on real or synthetic data by end of day.

**Day 2 (Aug 26 — today).** Stages 3–5: boundary enrichment test, population-specific classification, and the size-adjusted per-block SV-type enrichment scan. This is the statistical core of the project — get p-values/FDR-corrected results on the real dataset today, even if plots/report are rough.

**Day 3 (Aug 27).** Stages 6–8 in parallel across the team (population-cluster correlation, SV-based clustering + ARI, gene overlay for Maria/Alistair's duplication/inversion question). Start Stage 9 integration in the afternoon so there's a working summary table before the final day.

**Day 4 / buffer (Aug 28).** Finish Stage 9 (report + plots, stretch web view), run the end-to-end integration test below, prepare the demo/writeup. Treat this as buffer, not new-development time.

## Testing steps

- **Per-stage sanity checks on a small slice:** run Stages 0–9 on one chromosome (e.g., chr22) or a synthetic haploblock+SV set before running genome-wide; confirm intersection counts (Stage 2), classification counts (Stage 4), and enrichment p-value distributions (Stage 5) look sane (e.g., p-values roughly uniform under a shuffled-label negative control).
- **Negative controls:** re-run Stage 3's permutation test and Stage 5's regression on label-shuffled data — both should show no significant enrichment; a positive result on shuffled data indicates a pipeline bug, not a biological signal.
- **End-to-end integration run:** execute Stages 0→9 on the full nstd152 + haploblocks.org dataset for one chromosome, confirm the pipeline completes without manual intervention and the final summary table has one row per haploblock with all expected columns populated.
- **Validation against known structure:** in Stage 7, confirm the SV-based PCA/UMAP recovers *at least* the coarse 1KG superpopulation split (a well-established ground truth) before trusting its agreement with the haploblocks.org clusters specifically — if it can't recover superpopulations, the SV genotype matrix itself is broken, not the comparison.
- **Reproducibility check:** re-run the full pipeline twice with the same config (fixed seeds) and confirm identical output; re-run once with a different seed and confirm permutation/UMAP results are stable within expected tolerance.

## Non-pip dependencies

- `bcftools` / `samtools` (CLI) — VCF filtering/normalization in Stage 1.
- R is not required; the plan above is Python-only. If a team member prefers R for Stage 5/6 stats (e.g. `MASS::glm.nb`), that's a drop-in alternative to `statsmodels`, not a pipeline dependency.

## Prototype: Stage 0 (data ingestion)

`pipeline/stage0_ingest.py` is implemented and working — it's the current state of the prototype. **By default it fetches real data**: dbVar's nstd152 call+region VCFs and a (capped) slice of data.haploblocks.org's per-block cluster-hash TSVs. Any source that fails to download or parse falls back to synthetic data for just that input, so later stages are never blocked.

**Setup.** The `pyEnv_SVhack2026` virtualenv (`~/pyenvs/pyEnv_SVhack2026`) already has everything in `requirements.txt` installed. Elsewhere: `pip install -r requirements.txt` (plus `bcftools`/`samtools` for later stages).

**Run command** (no arguments needed — fetches real data by default, ~10-15s):
```
~/pyenvs/pyEnv_SVhack2026/bin/python pipeline/stage0_ingest.py
```
This writes `example_data/sv_calls.tsv`, `example_data/haploblocks.tsv`, `example_data/sample_metadata.tsv`, and `example_data/config.yaml`, plus the four raw dbVar files (`nstd152.GRCh38.variant_{call,region}.vcf.gz[.tbi]`). Useful flags:
- `--haploblock-chroms chr21,chr22` (default `chr21` only) and `--haploblock-max-blocks-per-chrom 50` (default) control how much of data.haploblocks.org's tens-of-thousands-of-files tree gets pulled — raise these, or pass `chroms all`, for a fuller run; expect it to be slow (one HTTP request per block file).
- `--skip-dbvar-download` / `--skip-haploblock-download` force synthetic data for that input, for fast offline iteration.
- `--sv-source` / `--haploblock-source` / `--panel-source` point at a single already-prepared local file or URL instead (bypasses the auto-fetchers above).

**Expected output** — log lines naming which of the three inputs are real vs. synthetic and key data-quality warnings, e.g.:
```
INFO: Using real haploblock table (50 blocks) from https://data.haploblocks.org/haploblock_hashes/1000G
INFO: dbVar call file: 102153 usable (sample, call) observation(s) across 9 sample(s)
INFO: Using real SV calls (50601 records) from ... (nstd152)
WARNING: No --panel-source given: replaced the placeholder synthetic sample panel with 9 sample(s)
  inferred from the real dbVar SV calls (hardcoded HGSVC trio lookup). This covers only
  ['AFR', 'AMR', 'EAS'] -- population-specific classification cannot be evaluated for any
  superpopulation not listed.
INFO: Wrote 50601 SVs, 50 haploblocks, 9 samples to .../example_data
```
`example_data/config.yaml` is the contract every later stage reads (genome build, AF/boundary thresholds, seeds, and the three table paths) — inspect it to confirm it parses and has all five top-level keys (`genome_build`, `data_sources`, `thresholds`, `seeds`, `paths`).

**Important data-availability caveat, not a bug:** nstd152's `variant_call.vcf.gz` currently covers exactly **9 samples from 3 HGSVC trios** (CHS: HG00512/13/14 → EAS, PUR: HG00731/32/33 → AMR, YRI: NA19238/39/40 → AFR) — no EUR or SAS samples at all. Without a real `--panel-source`, Stage 0 infers population labels for just these 9 via a hardcoded lookup (`HGSVC_TRIO_SUPERPOPULATIONS` in `stage0_ingest.py`). H3/H5/H6's population-level analyses will be severely underpowered on this alone (n=3 per population, 3 of 5 superpopulations represented) — a real IGSR sample panel widens the population *labels* available but does not add more dbVar samples, since nstd152 itself only has calls for these 9 individuals.

**Known limitations, deliberately out of scope for this prototype:** no retry/resume on failed downloads; per-block haploblock TSVs are parsed into summary stats (`n_clusters`, `hash_length`) but not persisted individually, since a genome-wide run touches tens of thousands of them; `cluster_diff_score` (an actual population-differentiation metric) is left as `NaN` for real haploblock data — computing it needs population labels, which belongs in a later stage (Stage 6), not ingestion; genome build is taken from explicit `--*-genome-build`/`--dbvar-build` flags rather than auto-detected from file headers.

`sv_calls.tsv` also carries `imprecise` (bool) and `length` (float, nullable) columns — dbVar never sets VCF QUAL/FILTER (always `"."` in nstd152), so `imprecise` is the closest thing this source has to a confidence flag, and INS records need `length` because their `end - start` is always ~1bp (an insertion's reference position, not its size); see the "Stage 1" comment block at the top of `stage0_ingest.py` for the full column contract.

## Prototype: Stage 1 (QC & normalization)

`pipeline/stage1_qc.py` is implemented and working. It reads Stage 0's `config.yaml`, applies a confidence filter (drops `imprecise` calls — dbVar's only real stand-in for FILTER=PASS), a size filter (`min_sv_length`/`max_sv_length` from config, INS rows with unresolvable length are exempted rather than dropped), coordinate sanity-checks (`start <= end`), and dedup, then validates the haploblock table is sorted and non-overlapping per chromosome — raising with the specific offending rows listed, not silently continuing, if it isn't.

**Dedup note, learned from the real data:** deduping on `(chrom, start, end, sv_type)` alone is wrong — 2439/50601 real rows collide that way, almost all INS, because dbVar gives every insertion the same ~1bp reference interval regardless of its actual size, so two textually "identical" rows are usually two distinct real insertions at the same locus. Adding `length` to the key resolves nearly all of that; the ~600 groups that still collide even with `length` included turn out to be the same event reported twice at different confidence (one precise, one imprecise) — ties are broken by keeping the precise record.

**Run command** (chains off Stage 0's output):
```
~/pyenvs/pyEnv_SVhack2026/bin/python pipeline/stage1_qc.py --config example_data/config.yaml --out-dir stage1_output
```
Writes `stage1_output/sv_calls.tsv`, `haploblocks.tsv`, `sample_metadata.tsv` (copied through), `qc_report.json` (before/after counts per filter step), and its own `config.yaml` (same schema as Stage 0's, `paths` repointed at these cleaned files) so Stage 2 can chain via `--config` the same way. `--min-sv-length`/`--max-sv-length`/`--keep-imprecise` override Stage 0's config values without re-running ingestion.

**Expected output** on the real chr21 + nstd152 data pulled by Stage 0's default run:
```
INFO: Confidence filter: dropped 14198 imprecise row(s) (drop_imprecise=True)
INFO: Size filter [50, 5000000]bp: dropped 10 row(s), exempted 2 with unknown length
INFO: Dedup: removed 30 exact-duplicate row(s)
INFO: Haploblock validation passed: 50 block(s), sorted and non-overlapping per chromosome
INFO: Wrote 36363/50601 SVs and 50 haploblocks to .../stage1_output
```

**Known limitations, deliberately out of scope:** true left-alignment (VCF-style indel normalization, shifting a breakpoint to the leftmost coordinate equivalent under the reference sequence) needs the reference genome; `--reference-fasta` is accepted but is an explicit no-op (logged loudly), not a best-guess implementation, since there's no reference data in this prototype to validate one against — `bcftools norm -f` is the standard tool for this if/when a reference is added.

## Prototype: Stage 2 (SV × haploblock intersection)

`pipeline/stage2_intersect.py` is implemented and working — **pure Python/numpy, no bedtools/pybedtools**. Haploblocks are non-overlapping within a chromosome (Stage 1 enforces this, and Stage 2 re-validates it defensively since it's independently runnable), so per chromosome the intersection is a sorted 1D interval search via `numpy.searchsorted` rather than a general-purpose interval-join library — fast (36,363 real SVs against 50 real haploblocks classify in well under a second) and dependency-free.

Each SV gets `position_class` (`within_block` / `boundary_crossing` / `outside_block`) and `haploblock_id` (comma-joined if it touches more than one block). A block only counts as fully "within" when the SV is ≥ N bp (`boundary_distance_bp`) from *both* its edges; real data.haploblocks.org blocks are contiguous (no gaps between neighbors, confirmed on the chr21 data pulled by Stage 0), so being near one block's edge usually means being near its neighbor's edge too — both get listed. See the classification rules spelled out in the module docstring for the exact cases (including the pandas empty-string-becomes-NaN gotcha in the round-tripped TSV, documented there for whoever implements Stage 3/4).

**Run command** (chains off Stage 1's output):
```
~/pyenvs/pyEnv_SVhack2026/bin/python pipeline/stage2_intersect.py --config stage1_output/config.yaml --out-dir stage2_output
```

**Expected output** on the real chr21 + nstd152 data:
```
INFO: Position classification (N=5000bp): within_block=19 boundary_crossing=4 outside_block=36340
INFO: 3 SV(s) matched more than one haploblock (span a shared edge)
INFO: Wrote 36363 annotated SVs and 50 haploblocks to .../stage2_output
```
(Most SVs are `outside_block` because Stage 0's default fetch only pulls 50 of chr21's ~6,099 haploblocks and every other chromosome has none at all — raise `--haploblock-max-blocks-per-chrom`/`--haploblock-chroms` in Stage 0 for fuller coverage.)

## Tests

`tests/test_stage0_ingest.py`, `tests/test_stage1_qc.py`, and `tests/test_stage2_intersect.py` run all three stages end-to-end on synthetic (offline) data — plus hand-crafted small-table cases for Stage 2's exact classification rules — and check the config/table contracts above. Run with:
```
~/pyenvs/pyEnv_SVhack2026/bin/python -m pytest tests/ -v
```
