# Stage 7: hash representation audit

## Purpose

Stage 7 directly evaluates what the haploblock representation captures or misses. It combines an SV-centric carrier-concentration view with a cluster-centric purity view, then records an interpretable category for each SV–haploblock pair.

## Inputs

The stage uses `sv_genotypes`, `sv_block_summary`, `sv_to_clusters`, `cluster_memberships`, and `sample_metadata` from the Stage 6 configuration. If Stage 4 classification is available, it also adds population context.

## Methods

### Carrier concentration

Stage 1's passing SV-to-cluster associations supply `expected_alt_haplotypes` as carrier evidence. The stage sums this evidence by SV–block–cluster. The leading cluster's share is the fraction assigned to the largest cluster. The effective carrier-cluster count is the inverse concentration of that evidence: it is near one if one cluster dominates and grows as evidence is distributed across clusters.

Both all-evidence and `standard`-evidence versions are reported. The standard version excludes low-evidence Stage 1 associations and is used for representation categories.

### Within-cluster purity

For every sufficiently represented cluster, Stage 7 counts called samples carrying at least one cluster haplotype and separates them into SV carriers and non-carriers. It reports the larger of those two proportions as `cluster_purity`; it reports their balance as `mixed_balance`. A mixed-cluster candidate must meet the configured minimum number of both carriers and non-carriers.

It also assesses complete local diplotypes. Because VCF phase is not linked to the haploblock `hap0`/`hap1` labels, a mixed complete diplotype is stronger evidence that the hash does not resolve the SV than a mixed single cluster alone.

### Representation categories

| `representation_pattern` | Rule |
|---|---|
| `hash_tag_candidate` | Exactly one standard-evidence carrier cluster, high carrier rate in that cluster, and no supported mixed diplotype. |
| `multi_cluster_sv_candidate` | More than one standard-evidence carrier cluster and no supported mixed diplotype. |
| `hash_subdivision_candidate` | One or no standard-evidence carrier cluster, but at least one supported mixed diplotype. |
| `multi_cluster_and_subdivision_candidate` | Multiple standard-evidence clusters and a supported mixed diplotype. |
| `insufficient_or_partial_evidence` | None of the above. |

The population-context flag `population_enriched_on_shared_cluster_candidate` identifies a Stage 4 population-specific SV whose top standard cluster occurs in at least the configured number of populations.

## Outputs

| File | Row grain | Main fields |
|---|---|---|
| `sv_carrier_cluster_summary.tsv` | SV–block | Carrier cluster count, leading cluster, leading share, and effective cluster count. |
| `sv_cluster_purity.tsv` | SV–block–cluster | Called samples, carrier/non-carrier counts, carrier rate, purity, balance, and mixed-count flag. |
| `sv_hash_representation.tsv` | SV–block | Consolidated representation and population-context patterns. |
| `sv_haploblock_information.tsv` | SV–block | Secondary information-gain metrics. |
| `haploblock_information_summary.tsv` | Block | Summary of secondary information-gain metrics. |
| `sv_pca_coordinates.tsv`, `sv_pca_variance.tsv`, `sv_pca.png` | Sample/component | Optional genotype-structure quality-control outputs. |

## Running and caveats

```bash
python claude-first-prototype/pipeline/stage7_information_gain.py \
  --config claude-first-prototype/stage6_output/config.yaml
```

Key defaults are four samples per cluster/diplotype, three carriers and three non-carriers for a mixed candidate, 0.9 purity for a tag candidate, and three populations for a cosmopolitan top cluster. `--skip-information-gain` and `--skip-pca` omit secondary diagnostics. PCA is quality control for cohort structure; it is not evidence that haploblocks represent SVs well. The categories are candidates requiring adequate support and replication, not final biological classifications.
