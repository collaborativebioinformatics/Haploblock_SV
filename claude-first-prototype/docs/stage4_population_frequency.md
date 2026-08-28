# Stage 4: population allele-frequency classification

## Purpose

Stage 4 provides population context independently of haploblock clusters. It prevents a simple co-occurrence between cluster and SV from being interpreted as a local relationship when both are common in the same population.

## Inputs and method

The stage reads `sv_genotypes` and `sample_metadata` from the Stage 1 configuration. For every SV and population, it counts called alternate alleles and divides by all called alleles to obtain alternate-allele frequency (AF). Missing genotypes do not contribute to either count.

The overall class uses only populations with at least `--min-samples-per-population` called samples. Defaults are a presence threshold of 0.05, an absence threshold of 0.01, and two called samples per population.

| Class | Rule |
|---|---|
| `common` | AF meets the presence threshold in at least two populations. |
| `population_specific` | AF meets it in exactly one population and is below the absence threshold in all other populations with data. |
| `other` | Insufficient data, rare/absent, or an intermediate pattern not captured by the two classes. |

## Outputs

| File | Row grain | Important fields |
|---|---|---|
| `sv_af_classification.tsv` | SV–population | Fixed SV metadata, `population`, `n_samples`, `n_called`, `called_alleles`, `pop_has_data`, `af`. |
| `sv_classification.tsv` | SV | Fixed SV metadata, `sv_class`, `specific_to_population`, `other_reason`. |
| `config.yaml` | Run | Carried-forward Stage 1 paths plus the frequency settings and both output paths. |

`n_samples` is the number assigned to the population; `n_called` is the subset with a called genotype for that SV. `called_alleles` is normally twice the called-sample count for diploid calls. `af` is the proportion of those alleles carrying the alternate allele.

## Running and caveat

```bash
python claude-first-prototype/pipeline/stage4_classify_af.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

Class labels are descriptive and depend strongly on sample size and thresholds. They are not tests of population differentiation.
