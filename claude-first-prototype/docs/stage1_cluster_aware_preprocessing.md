# Stage 1: cluster-aware preprocessing

## Purpose

Stage 1 is the core data-preparation and association stage. It reads the cohort SV VCF, aligns its samples to haploblock membership data, and asks which haploblock clusters support each SV. It writes the standardized tables consumed by all later analyses.

The result is an **SV–haploblock–cluster association**, not simply an overlap with a fixed genomic region. Each overlapping SV is evaluated independently against every cluster in that block. A record can have no passing clusters, several passing clusters, or passing clusters in several blocks.

## Inputs

| Input | Default or required content | Use |
|---|---|---|
| Cohort VCF | `input/1kgp_ont_cohort.postfilter.full.vcf.gz`; GRCh38 | Source SV coordinates, metadata, and per-sample genotypes. |
| Haploblock cluster files | Downloaded from the 1000G haploblock-hash endpoint or supplied with `--cluster-root` | Defines blocks and sample-haplotype cluster memberships. |
| Sample metadata | Downloaded 1000 Genomes ONT metadata or `--sample-metadata` | Adds population labels for Stages 4, 6, and 7. |
| GTF annotation | Ensembl GRCh38 release 115 or `--gtf` | Registered for Stage 8 annotation. |

The VCF and cluster files must describe compatible sample identities. Stage 1 changes `GM<number>` sample IDs to `NA<number>` and writes both identities to `samples.tsv`. Metadata can use either form; the output uses canonical IDs.

## Coordinate and genotype conventions

All written SV and block coordinates are 0-based, half-open intervals: the start is included and the end is excluded. An SV overlaps a block when `sv_start < block_end` and `sv_end > block_start`.

`sv_record_id` is a new unique identifier in the form `chrN_record_k`, scoped to the run. Use it for joins and counting. `sv_id` is retained from the VCF and can be repeated or shortened; it is not a primary key. A missing VCF ID is replaced with a readable source label, and overly long source labels are shortened with a stable hash suffix.

Genotypes containing `.` are treated as missing rather than reference. Stages consuming dosage interpret alternate allele dosage as 0, 1, or 2 for called diploid genotypes.

## Association method

For each SV that overlaps a block, Stage 1 uses the samples whose haplotypes are represented in the cluster file. Homozygous alternate and reference genotypes give unambiguous allele evidence. For a heterozygote, the VCF phase is not assumed to correspond to the cluster file's `hap0`/`hap1` labels. The method instead uses the evidence from the two relevant clusters to estimate the posterior probability that the alternate allele belongs to either cluster.

An expectation-maximization (EM) calculation estimates the probability that a haplotype in each cluster carries the SV. A cluster is written to `sv_to_clusters` when the estimate meets `--association-threshold` and there is enough callable evidence. A heterozygote contributes an allele to one cluster only if its posterior assignment reaches `--posterior-threshold`; otherwise it remains ambiguous.

The default thresholds are:

| Parameter | Default | Meaning |
|---|---:|---|
| `--association-threshold` | 0.75 | Minimum estimated probability of SV carriage for a cluster to be reported. |
| `--posterior-threshold` | 0.75 | Minimum confidence needed to assign a heterozygous allele to one cluster. |
| Required callable haplotypes | min(3, represented haplotypes) | Adaptive minimum evidence for the cluster estimate. |
| `--max-iterations` | 25 | Maximum EM updates. |
| `--tolerance` | 1e-5 | Change below which the EM estimate is considered converged. |

Cluster `hap0` and `hap1` labels are technical labels, not maternal and paternal labels.

## Outputs and data contract

| File | Row grain | Key fields | Downstream role |
|---|---|---|---|
| `sv_genotypes.<chrom>.tsv` | One VCF record | `sv_record_id`; fixed SV metadata; one raw `GT` column per canonical sample | Genotype contract for Stages 4, 6, and 7. |
| `samples.tsv` | One VCF sample | `sample_id`, `original_sample_id` | Maps canonical IDs back to the source VCF. |
| `sample_metadata.tsv` | One represented sample | `sample_id`, `population`, optional `superpopulation` | Population contract. Rows absent from the VCF are removed and order follows the VCF. |
| `haploblocks.<chrom>.tsv` | One block | `haploblock_id`, `chrom`, `start`, `end` | Coordinates and denominator for spatial/count analyses. |
| `cluster_memberships.<chrom>.tsv` | One sample haplotype in a block cluster | `haploblock_id`, `sample_id`, `haplotype`, `cluster_id` | Complete cluster universe for Stages 6 and 7. |
| `sv_block_summary.<chrom>.tsv` | One overlapping SV–block pair | `sv_record_id`, `haploblock_id` | Candidate universe; never restricts to passing associations. |
| `sv_to_clusters.<chrom>.tsv` | One passing SV–block–cluster association | `sv_record_id`, `haploblock_id`, `cluster_id` | Supported Stage 1 associations and carrier-concentration evidence. |
| `config.yaml` | One run | Paths, source information, and thresholds | Input contract for later scripts. |

### Fixed SV metadata columns

All principal SV tables begin with these columns:

| Field | Meaning |
|---|---|
| `sv_record_id` | Unique Stage 1 record key. |
| `sv_id` | Source VCF ID or a derived readable label; not guaranteed unique. |
| `chrom`, `start`, `end` | GRCh38 SV interval in 0-based, half-open coordinates. |
| `sv_type`, `length` | SV type and non-negative length derived from `SVLEN` or interval length. |
| `filter`, `imprecise` | VCF `FILTER` value and whether the INFO field contains `IMPRECISE`. |

### `sv_block_summary` fields

| Field | Meaning |
|---|---|
| `block_start`, `block_end` | Coordinates of the overlapping haploblock. |
| `overlaps_multiple_blocks` | True when the SV overlaps more than one haploblock. |
| `callable_samples`, `heterozygous_samples` | SV evidence available among represented samples. |
| `clusters_evaluated`, `associated_clusters`, `associated_cluster_ids` | Number assessed, number passing, and IDs of passing clusters. |
| `association_class` | Summary of the Stage 1 association outcome for the SV–block pair. |
| `model_converged`, `model_iterations` | EM convergence context. |

### `sv_to_clusters` fields

| Field | Meaning |
|---|---|
| `cluster_id` | Local haploblock-cluster identifier. |
| `cluster_haplotypes_total`, `cluster_haplotypes_in_vcf` | Haplotype count in the source cluster file and the subset represented in the VCF. |
| `callable_haplotypes`, `required_callable_haplotypes`, `call_rate` | Amount and fraction of usable VCF evidence. |
| `expected_alt_haplotypes` | Alternate-bearing haplotype evidence, including probabilistic assignments for heterozygotes. |
| `sv_probability` | EM estimate that a haplotype in this cluster carries the SV. |
| `ci95_low`, `ci95_high` | Approximate uncertainty interval for the estimate. |
| `evidence_tier` | `low` below three callable haplotypes; otherwise `standard`. |
| `model_converged`, `model_iterations` | EM convergence context shared by the SV–block estimate. |

An SV absent from `sv_to_clusters` has no association that passed the configured rule. It is not evidence that the SV is absent from all clusters.

## Configuration

`config.yaml` records the genome build, input sources, thresholds, run settings, and absolute paths to every Stage 1 output. Downstream scripts should use paths from this configuration rather than infer filenames. Stage 4 and later carry it forward and add their own output paths and settings.

## Running Stage 1

```bash
python claude-first-prototype/pipeline/stage1_cluster_aware.py \
  --vcf input/1kgp_ont_cohort.postfilter.full.vcf.gz \
  --chroms chr21 \
  --out-dir claude-first-prototype/stage1_output
```

Useful options are `--cluster-root` for local cluster files, `--sample-metadata` and `--gtf` for local reference inputs, and `--threads` for independent chromosome or block work. Downloaded cluster files are cached under `_intermediate/clusters/` by source URL. The VCF is indexed with Tabix when possible; otherwise it is streamed.

## Caveats

Stage 1 preserves input calls and reports their quality. It does not resolve whether distinct records in a repetitive locus are the same biological allele—that belongs to Stage 0. A low-evidence association can be useful for inspection but should not be used as definitive carrier evidence. An association in several blocks reflects a spanning SV being evaluated in each overlapping local context, not duplicated biological events.
