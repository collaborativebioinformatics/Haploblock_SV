# Stage 6: population-conditioned SV–cluster association

## Purpose

Stage 6 tests whether a local haploblock cluster predicts SV carriage beyond population membership. It is the primary test for a cluster that might be a portable proxy for direct SV genotyping.

## Inputs

The stage reads `sv_genotypes`, `sv_block_summary`, `cluster_memberships`, and `sample_metadata` from a carried-forward configuration (normally Stage 4). It tests every eligible cluster in every block overlapped by each SV, including clusters that did not pass Stage 1's association threshold.

## Method

SV dosage is 0, 1, or 2 alternate alleles per sample. Cluster dosage is the number of a sample's haplotypes assigned to that cluster (0, 1, or 2). For each population separately, the stage subtracts the population mean from both quantities, then calculates their correlation. This `population_adjusted_r` asks whether samples differ together **within** their population.

For an empirical null, SV dosages are shuffled only among samples from the same population. This keeps the observed population frequency pattern intact. Within an SV–block pair, the permutation procedure uses the largest cluster statistic to account for the multiple clusters assessed in that block. Benjamini–Hochberg adjustment is then applied across SV–block pairs.

The run triages all pairs with a small number of permutations, screens promising pairs more deeply, and independently refines the most promising results. Results below `--min-abs-r` do not advance to expensive stages but remain in the association output.

## Outputs

| File | Row grain | Key fields |
|---|---|---|
| `sv_cluster_associations.tsv` | Tested SV–block–cluster | Fixed SV metadata, `haploblock_id`, `cluster_id`, counts/rates, effect, probability values, and cross-population consistency. |
| `sv_cluster_summary.tsv` | SV–block | Fixed SV metadata, best cluster IDs, best effect and q-value, and `association_pattern`. |
| `config.yaml` | Run | Carried-forward paths, association settings, and output paths. |

### Association fields

| Field | Meaning |
|---|---|
| `n_called` | Samples with a called SV dosage in the tested block. |
| `n_cluster_carriers`, `n_samples_with_cluster` | Cluster-haplotype and sample support. |
| `n_sv_carriers_with_cluster`, `n_sv_noncarriers_with_cluster` | Direct carrier/non-carrier counts among samples carrying the cluster. |
| `carrier_rate_with_cluster`, `carrier_rate_without_cluster`, `carrier_rate_difference` | Unadjusted, easy-to-inspect prevalence difference. |
| `association_direction` | `carrier_enriched` or `carrier_depleted`. A depleted cluster is exclusion evidence, not a carrier tag. |
| `population_adjusted_r` | Direction and strength after population averages are removed. |
| `p_value`, `q_value` | Permutation probability and across-pair adjusted probability. |
| `permutations_used` | Permutations used for the retained result, not the total spent on discarded triage work. |
| `informative_populations`, `directional_consistency` | Number of populations with usable within-population comparisons and fraction matching the overall direction. |

`association_pattern` is `cross_population_consistent_tag_candidate` when an enriched result is significant, sufficiently large, and directionally consistent across at least two informative populations. Other labels distinguish population-dependent association, a detected association, an exclusion signal, and no detected signal.

## Running and caveats

```bash
python claude-first-prototype/pipeline/stage6_cluster_association.py \
  --config claude-first-prototype/stage4_output/config.yaml
```

Defaults include 10,000 screening permutations, 1,000,000 refinement permutations, a minimum cluster size of six haplotypes, a minimum effect magnitude of 0.3, and q-value threshold 0.05. Directional consistency is an in-sample screen, not independent evidence that the tag will transfer to another cohort. Batch effects and sparse populations can still produce misleading candidates.
