# Stage 4: population allele-frequency classification

## Purpose

Stage 4 provides population context independently of haploblock clusters. It prevents a simple co-occurrence between cluster and SV from being interpreted as a local relationship when both are common in the same population.

Stage 4 uses population labels independently of haploblock cluster labels.
This helps distinguish a local SV-haploblock association from a pattern that
may be explained simply by population allele-frequency differences.

## Inputs and method

The stage reads `sv_genotypes` and `sample_metadata` from the Stage 1 configuration. For every SV and population, it counts called alternate alleles and divides by all called alleles to obtain alternate-allele frequency (AF). Missing genotypes do not contribute to either count.

The Stage 4 `population` field represents the 1000 Genomes subpopulation
derived from the metadata `SubPopulation` field, for example `YRI`, `PEL`,
`GBR`, or `ACB`. The metadata may also contain a `superpopulation` field,
such as `AFR`, `EUR`, or `AMR`, but Stage 4 currently groups samples using
`population`, not `superpopulation`.

The overall class uses only populations with at least `--min-samples-per-population` called samples. Defaults are a presence threshold of 0.05, an absence threshold of 0.01, and two called samples per population.

| Class | Rule |
|---|---|
| `common` | AF meets the presence threshold in at least two populations. |
| `population_specific` | AF meets it in exactly one population and is below the absence threshold in all other populations with data. |
| `other` | Insufficient data, rare/absent, or an intermediate pattern not captured by the two classes. |

## Outputs

Stage 4 writes population-level allele-frequency results, one-row-per-SV
classifications, and haploblock-level summaries.

| File | Row grain | Important fields |
|---|---|---|
| `sv_af_classification.tsv` | One SV-population pair | Fixed SV metadata, `population`, `n_samples`, `n_called`, `called_alleles`, `pop_has_data`, `af`. |
| `sv_classification.tsv` | One row per SV | Fixed SV metadata, `sv_class`, `specific_to_population`, `other_reason`. |
| `sv_classification_haploblocks.tsv` | One SV-haploblock pair | Stage 4 classification fields plus `haploblock_id`, `block_start`, and `block_end`. |
| `stage4_summary.tsv` | One row per run | Total unique SVs, common SVs, population-specific SVs, other SVs, and haploblock overlap counts. |
| `population_specific_summary.tsv` | One row per subpopulation | Number of population-specific SVs and number of haploblocks containing them. |
| `haploblock_population_specific_summary.tsv` | One row per haploblock | `total_svs`, `common_svs`, `population_specific_svs`, populations with population-specific SVs, and block coordinates. |
| `config.yaml` | One row per run | Carried-forward Stage 1 paths, Stage 4 thresholds, and paths to the Stage 4 output files. |
| `stage4_plots/<chrom>/` | One directory per chromosome | PNG figures and TSV files summarizing SV classes, population-specific SVs, and haploblock position. |

`n_samples` is the number of samples assigned to a subpopulation. `n_called`
is the number of samples with a non-missing genotype for that SV.
`called_alleles` is the number of called alleles used as the AF denominator.
For diploid genotypes, this is normally twice the number of called samples.
`af` is the alternate-allele frequency calculated as alternate alleles divided
by called alleles.

`sv_classification_haploblocks.tsv` is created by joining
`sv_classification.tsv` to Stage 1's `sv_block_summary` using
`sv_record_id`. A populated `haploblock_id` identifies an SV-haploblock
overlap. An empty value means that the SV does not overlap a haploblock.

One SV can overlap more than one haploblock, so
`sv_classification_haploblocks.tsv` can contain more rows than
`sv_classification.tsv`. The haploblock summary counts unique
`sv_record_id` values within each haploblock to prevent duplicate
SV-haploblock records from being counted more than once.

Use `sv_record_id` as the join key. It uniquely identifies one input VCF
record, whereas `sv_id` is a source label that may be repeated by the variant
caller.

The primary output for identifying haploblocks containing
population-specific SVs is `haploblock_population_specific_summary.tsv`.
Rows are ranked by the number of population-specific SVs in each haploblock.
The detailed `sv_classification_haploblocks.tsv` file should be used to
inspect the individual SVs, their assigned subpopulations, and their
haploblock coordinates.

### Plot outputs

Plots are generated only when the `--plots` option is supplied. Each
chromosome listed in the Stage 1 configuration receives its own subdirectory,
so a run containing `chr21` and `chr22` produces separate plot directories:


```text
stage4_plots/
├── chr21/
└── chr22/
```
The current plot outputs are:

```text
stage4_plots/<chrom>/
├── haploblock_sv_classification.png
├── population_specific_by_haploblock.png
├── population_specific_fraction.png
├── haploblock_plot_data.tsv
└── population_haploblock_plot_data.tsv
```



## Running and caveat

```bash
python claude-first-prototype/pipeline/stage4_classify_af.py \
  --config claude-first-prototype/stage1_output/config.yaml
```

Class labels are descriptive and depend strongly on sample size and thresholds. They are not tests of population differentiation.
