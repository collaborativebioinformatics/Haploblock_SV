# Haploblock_SV: structural variants within haploblocks

Haploblock_SV asks a practical biological question: **when can a local, SNV-derived haploblock cluster stand in for direct structural-variant (SV) genotyping, and when does the SV add information that the cluster misses?**

Haploblocks, defined by genomic hashing in the [haploblocks.org](https://haploblocks.org) project, group similar local haplotypes. Structural variants can be difficult to genotype at scale, especially in cohorts without long-read sequencing. If an SV is consistently carried by one haploblock cluster, the cluster may be a useful proxy. If carriers and non-carriers are mixed within the same cluster, direct SV genotyping adds resolution.

The pipeline begins with a merged cohort SV VCF and produces interpretable, locus-level tables. It uses population labels to avoid mistaking broad ancestry differences for a local haploblock–SV relationship.

For inputs, calculations, parameters, and complete output contracts, see the [stage-by-stage methods documentation](docs/README.md).

## What the results mean

Two complementary views describe whether a haploblock represents an SV well:

| View | Question | A result that supports the haploblock | A result that suggests the SV adds information |
|---|---|---|---|
| SV-centric | Where are carriers of this SV found? | Carrier evidence is concentrated in one cluster. | Carriers occur in several clusters. |
| Cluster-centric | Are members of this cluster alike for SV carriage? | Nearly all callable members are carriers, or nearly all are non-carriers. | The same sufficiently represented cluster contains both carriers and non-carriers. |

These are association results, not proof of causality or functional effect. A gene overlap identifies a locus worth examining; it does not establish that the SV changes that gene's function.

## Pipeline at a glance

Stages are numbered to match the scripts. Stage 3 has no user-facing analysis.

| Stage | Biological question | Main result |
|---|---|---|
| 0 | Are equivalent SV calls represented as one comparable cohort allele? | A merged, reconciled cohort VCF (planned). |
| 1 | Which haploblock clusters are compatible with each SV? | The reusable SV, haploblock, and cluster tables used by downstream stages. |
| 2 (optional) | Does an SV sit within a block, cross a boundary, or lie near one? | Descriptive boundary classes. |
| 4 | Is an SV common across populations or concentrated in one? | Per-population allele frequencies and an SV classification. |
| 5 (optional) | Does a block contain more of one SV type than expected for its length? | A screening table of block-by-SV-type enrichments. |
| 6 | Does a cluster predict SV carriage after accounting for population? | Candidate cluster tags, population-dependent signals, and exclusion signals. |
| 7 | Does the haploblock capture, or miss, SV information? | Carrier concentration, within-cluster purity, and representation categories. |
| 8 | Which genes or exons lie at credible captured or missed SV candidates? | Type-aware gene and breakpoint annotations. |
| 9 | What is the biological summary across stages? | A concise Markdown/HTML report and figures. |

## Data scope and important caveats

The current prototype starts with `input/1kgp_ont_cohort.postfilter.full.vcf.gz`, a previously merged 1000 Genomes ONT callset. It contains deletions and insertions because post-collapse kanpig regenotyping required resolved variant sequences; inversions and breakends are not present in this input. This is a limitation of the current callset, not of the biological question.

The VCF phase is not assumed to correspond to haploblock `hap0` and `hap1` labels. Consequently, an unphased heterozygote is assigned to a cluster only when the surrounding evidence supports that assignment; otherwise it remains ambiguous. Haploblock haplotype labels are not maternal and paternal labels.

Rare calls, missing genotypes, SV representation errors in repetitive sequence, sequencing batches, and uneven population sampling can all create apparent associations. Treat candidate tags and candidate within-cluster splits as prioritization results until they have adequate support and replication.

## Stage 0 — cohort SV merging

### Biological rationale

The same biological SV can be reported with slightly different breakpoints or representations by different samples and callers. Reconciling those records is essential before interpreting an allele on several haplotype backgrounds: otherwise an apparent multi-background allele may simply be a collection of non-equivalent calls.

### Outputs

## Stage 0: cohort SV merging

[`pipeline/stage_0_mergingvcf.py`](pipeline/stage_0_mergingvcf.py) accepts a
manifest of local single-sample VCFs, validates sample names and contig
compatibility, sorts and uniquely identifies each input, merges them with
`bcftools merge -m id`, and reconciles equivalent representations with
`truvari collapse`. See [`STAGE0_SUMMARY.md`](STAGE0_SUMMARY.md) for the command, 
output contract, and options.

| Output | Contents | Status |
|---|---|---|
| Merged cohort VCF | One cohort-level representation of comparable SV calls, with sample genotypes. | Planned; not implemented in this prototype. |

### Implementation notes


## Stage 1 — cluster-aware preprocessing

### Biological rationale

Stage 1 creates an auditable mapping between cohort SVs and local haploblock clusters. This is the shared evidence needed to decide whether a cluster tags an SV, whether an SV occurs on several backgrounds, or whether a cluster contains mixed carriers and non-carriers.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `sv_genotypes.<chrom>.tsv` | Input VCF record | Normalized SV metadata and raw genotypes; the source for population and representation analyses. |
| `samples.tsv` | Sample | Canonical sample ID and original VCF sample ID. |
| `sample_metadata.tsv` | Represented sample | Population metadata, ordered to the VCF. |
| `haploblocks.<chrom>.tsv` | Haploblock | Block coordinates and IDs. |
| `cluster_memberships.<chrom>.tsv` | Sample haplotype in a cluster | The complete cluster membership universe. |
| `sv_block_summary.<chrom>.tsv` | Overlapping SV–haploblock pair | Every SV overlapping a block; use for candidate-universe analyses. |
| `sv_to_clusters.<chrom>.tsv` | Passing SV–haploblock–cluster association | Clusters that pass Stage 1 support and probability thresholds. |
| `config.yaml` | Run | Registered paths and run settings for the later stages. |

### Key output fields

| Table | Fields | Meaning |
|---|---|---|
| `sv_genotypes` | `sv_record_id`, `sv_id`, `chrom`, `start`, `end`, `sv_type`, `length`, `filter`, `imprecise` | `sv_record_id` is the unique join key. `sv_id` is a source label and may not be unique. Coordinates are 0-based, half-open. Sample genotype columns follow the metadata columns. |
| `cluster_memberships` | `haploblock_id`, `sample_id`, `haplotype`, `cluster_id` | Identifies the two haplotypes assigned to a sample and their local cluster; use it to enumerate all clusters in a block. |
| `sv_block_summary` | `sv_record_id`, `haploblock_id` | A unique SV–block pair, including pairs with no passing Stage 1 cluster association. |
| `sv_to_clusters` | `cluster_id`, `cluster_haplotypes_in_vcf`, `callable_haplotypes`, `call_rate`, `expected_alt_haplotypes`, `sv_probability`, `ci95_low`, `ci95_high`, `evidence_tier` | A passing association. `sv_probability` estimates the chance that a haplotype in the cluster carries the SV. The interval shows uncertainty around that estimate; `low` evidence means fewer than three callable haplotypes. |
| `sample_metadata` | `sample_id`, `population`, optional `superpopulation` | Population is required for Stages 4 and 6. |

Absence from `sv_to_clusters` means no cluster passed the configured rule; it does **not** show that the SV is absent from every cluster.

### Implementation notes

Run `python claude-first-prototype/pipeline/stage1_cluster_aware.py`. Stage 1 reads the merged VCF, retrieves or reuses haploblock cluster membership, and can retrieve matching 1000 Genomes metadata and Ensembl GRCh38 annotation. Cluster files are cached under `_intermediate/clusters/`; use the generated `config.yaml` as the input to downstream scripts.

## Stage 2 — boundary classification (optional)

### Biological rationale

An SV that crosses a haploblock boundary may identify a locus where the fixed block boundaries are not a good spatial description of structural variation. This is a distinct spatial observation, not evidence that a cluster predicts the SV.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `boundary_svs.<chrom>.tsv` | SV record | Position relative to haploblock boundaries. |
| `boundary_qc.json` | Run and chromosome | Counts of position classes and near-boundary records. |

### Key output fields

| Field | Meaning |
|---|---|
| `position_class` | `outside_block` (no overlapping blocks), `within_block` (one), or `boundary_crossing` (two or more). |
| `near_boundary` | Whether the SV lies within the configured distance of a boundary; this is separate from crossing it. |
| `haploblock_id` | The overlapping block or blocks; empty for an SV outside all blocks. The distance threshold used for `near_boundary` is recorded in `boundary_qc.json`. |

### Implementation notes

Run `python claude-first-prototype/pipeline/stage2_intersect.py --config claude-first-prototype/stage1_output/config.yaml`. It reuses Stage 1 metadata rather than rescanning the VCF.

## Stage 4 — population allele-frequency classification

### Biological rationale

Population frequency provides essential context. A cluster and an SV may appear associated simply because both are common in the same population. This stage also identifies population-enriched SVs that occur on an otherwise widespread cluster, a pattern consistent with an event that the older SNV-defined cluster does not resolve.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `sv_af_classification.tsv` | SV–population | Population-specific allele-frequency evidence. |
| `sv_classification.tsv` | SV | Overall `common`, `population_specific`, or `other` class. |

### Key output fields

| Table | Fields | Meaning |
|---|---|---|
| `sv_af_classification` | `population`, `n_samples`, `n_called`, `called_alleles`, `pop_has_data`, `af` | `af` is the fraction of called alleles carrying the alternate SV. `pop_has_data` identifies populations meeting the minimum sample requirement. |
| `sv_classification` | `sv_class`, `specific_to_population`, `other_reason` | `common` meets the frequency threshold in at least two populations; `population_specific` meets it in one population and is below the absence threshold in the others; `other` is not cleanly classified. |

### Implementation notes

Run `python claude-first-prototype/pipeline/stage4_classify_af.py --config claude-first-prototype/stage1_output/config.yaml`. Default thresholds are 5% for presence, 1% for absence elsewhere, and two called samples per population; choose thresholds appropriate to the cohort size.

## Stage 5 — per-haploblock SV-type enrichment (optional)

### Biological rationale

An excess of one SV type in a block may point to local sequence architecture that promotes a particular mutational mechanism. It is a hypothesis-generating screen and does not test whether haploblock clusters represent SV alleles.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `sv_type_enrichment.tsv` | Haploblock–SV-type combination | Observed and length-adjusted expected counts, with a multiple-testing-adjusted flag. |

### Key output fields

| Field | Meaning |
|---|---|
| `observed_count` | Number of unique SV–block pairs of this type. |
| `expected_count` | Count expected if this SV type were distributed in proportion to block length. |
| `p_value` | How surprising the observed count is under that simple length-based model. |
| `q_value` | The p-value adjusted across all tested block/type combinations, controlling the expected proportion of false-positive flags. |
| `flagged` | Whether `q_value` is below the selected threshold. |

### Implementation notes

Run `python claude-first-prototype/pipeline/stage5_type_enrichment.py --config claude-first-prototype/stage1_output/config.yaml`. The model adjusts only for block length. Callability, SNP density, uneven variance in SV counts, and minimum-count choices should be considered before interpreting a flagged cell as a biological enrichment.

## Stage 6 — population-conditioned SV–cluster association

### Biological rationale

This stage tests whether a local cluster is associated with SV carriage **within populations**, rather than merely sharing a population distribution with the SV. A carrier-enriched cluster with a consistent direction in multiple populations is a candidate portable tag. A cluster depleted of carriers can provide useful exclusion information, but is not a carrier tag.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `sv_cluster_associations.tsv` | Tested SV–haploblock–cluster combination | Full association evidence for every eligible cluster. |
| `sv_cluster_summary.tsv` | SV–haploblock pair | Best cluster and an interpretable association category. |

### Key output fields

| Table | Fields | Meaning |
|---|---|---|
| `sv_cluster_associations` | `n_called`, `n_samples_with_cluster`, `n_sv_carriers_with_cluster`, `n_sv_noncarriers_with_cluster`, `carrier_rate_with_cluster`, `carrier_rate_without_cluster`, `carrier_rate_difference` | Direct counts and rates. A positive rate difference means carriers are more frequent among samples with the cluster. |
| `sv_cluster_associations` | `association_direction`, `population_adjusted_r`, `p_value`, `q_value`, `permutations_used`, `informative_populations`, `directional_consistency` | `population_adjusted_r` measures the association after removing population averages. `p_value` comes from shuffling SV genotypes within each population; `q_value` adjusts for the many SV–block tests. `directional_consistency` is the fraction of informative populations with the same direction—not independent validation. |
| `sv_cluster_summary` | `best_cluster_id`, `best_enriched_cluster_id`, `best_depleted_cluster_id`, `association_pattern` | `association_pattern` is `cross_population_consistent_tag_candidate`, `population_dependent_association`, `cluster_associated`, `cluster_exclusion_signal`, or `no_detected_cluster_signal`. |

### Implementation notes

Run `python claude-first-prototype/pipeline/stage6_cluster_association.py --config claude-first-prototype/stage4_output/config.yaml`. The procedure shuffles genotypes only within populations to retain the observed population frequency pattern. It first screens pairs with fewer permutations and refines promising ones; `permutations_used` is the count used for the retained result. Population sign consistency is an in-sample screen, not proof that a tag transfers to a new cohort.

## Stage 7 — hash representation audit

### Biological rationale

Stage 7 is the direct representation test. It measures both where carrier evidence lies across clusters and whether clusters themselves are internally consistent for carriage. It can identify a well-tagged SV, a candidate multi-cluster SV, or a count-supported candidate where the SV subdivides an existing hash group.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `sv_carrier_cluster_summary.tsv` | SV–haploblock pair | Concentration of supported carrier evidence across clusters. |
| `sv_cluster_purity.tsv` | SV–haploblock–cluster combination | Carrier and non-carrier composition of a cluster. |
| `sv_hash_representation.tsv` | SV–haploblock pair | Main representation category and population context. |
| `sv_haploblock_information.tsv` and `haploblock_information_summary.tsv` | SV–block and block | Secondary information-gain descriptions. |
| `sv_pca_coordinates.tsv`, `sv_pca_variance.tsv`, and PCA plot | Sample and component | Structure quality control only. |

### Key output fields

| Table | Fields | Meaning |
|---|---|---|
| `sv_carrier_cluster_summary` | `n_supported_carrier_clusters`, `top_supported_cluster_id`, `top_cluster_carrier_evidence_share`, `effective_carrier_cluster_count` | The top-cluster share is the proportion of supported carrier evidence in the leading cluster. An effective count near one means evidence is concentrated; larger values mean it is spread across clusters. `standard`-evidence versions of these fields exclude sparse associations. |
| `sv_cluster_purity` | `n_called_cluster_samples`, `n_sv_carriers`, `n_sv_noncarriers`, `carrier_rate_in_cluster`, `cluster_purity`, `mixed_balance`, `meets_mixed_count_threshold` | `cluster_purity` approaches 1 when a cluster is mostly one state (carrier or non-carrier). `mixed_balance` approaches 1 when the two states are evenly represented. The threshold flag requires enough of both states to call the cluster a candidate split. |
| `sv_hash_representation` | `representation_pattern`, `n_mixed_clusters_meeting_count_threshold`, `n_mixed_diplotypes_meeting_count_threshold`, `top_standard_cluster_carrier_rate` | The main biological interpretation. Mixed complete diplotypes are the stronger within-group evidence because the VCF phase is not tied to haploblock haplotype labels. |
| `sv_hash_representation` | `sv_class`, `specific_to_population`, `top_standard_cluster_population_count`, `top_standard_cluster_populations`, `population_context_pattern` | Describes the population context, including `population_enriched_on_shared_cluster_candidate` where applicable. |

### Implementation notes

Run `python claude-first-prototype/pipeline/stage7_information_gain.py --config claude-first-prototype/stage6_output/config.yaml --threads 4`. PCA and information gain are secondary diagnostics; PCA shows cohort structure and batch/ancestry patterns, not whether a haploblock adequately represents an SV. `--skip-information-gain` and `--skip-pca` omit those optional diagnostics while retaining the primary concentration and purity results.

## Stage 8 — consequence-aware candidate annotation

### Biological rationale

This stage gives the clearest captured and missed SV candidates gene, exon, and breakpoint context. It distinguishes consequences by SV type: for example, an inversion breakpoint can disrupt a gene even when genes contained within the inverted span are not directly interrupted.

### Outputs

| Output | One row per | Use |
|---|---|---|
| `annotated_sv_candidates.tsv` | Candidate SV–haploblock pair | Association, representation, call-quality, and type-aware annotation evidence kept in separate columns. |

### Key output fields

| Field | Meaning |
|---|---|
| `consequence` | Type-aware annotation category, such as breakpoint disruption, exon overlap, gene overlap, or a gene contained within an inversion. |
| `genes` | Overlapping or breakpoint-associated gene names, depending on `overlap_basis`. |
| `overlap_basis` | States whether the annotation comes from span overlap or a breakpoint. |
| `representation_pattern`, `population_context_pattern` | Stage 7 context explaining why this SV–block pair was selected. |
| `sv_class`, `specific_to_population` | Stage 4 population context. |
| `call_quality` | Combined VCF filter and imprecision label: `pass_precise`, `pass_imprecise`, `nonpass_precise`, or `nonpass_imprecise`. |
| `consequence_priority` | A technical sort key for the severity ordering of annotation categories; it is not a biological impact score. |

### Implementation notes

Run `python claude-first-prototype/pipeline/stage8_candidate_annotation.py --config claude-first-prototype/stage7_output/config.yaml`. It uses the Ensembl GRCh38 GTF registered by Stage 1 unless `--gtf` supplies another annotation. Do not collapse association strength, call quality, and functional annotation into one score: they are different kinds of evidence.

## Stage 9 — report and biological interpretation

### Biological rationale

The final report brings together the representation results and annotations so that users can distinguish robustly supported patterns from plausible biological hypotheses. It is intended to prioritize loci for follow-up, not replace read-level review or replication.

### Outputs

| Output | Contents |
|---|---|
| `report.md` and `report.html` | A factual summary of stages 5–8, with figures. |
| `report_facts.json` | Exact compact values used to build the report. |
| `figures/` | Enrichment, representation, and PCA summary plots when source tables are available. |
| `agent_metadata.json` | Whether an optional plain-language interpretation was added. |

### Key output fields

| Output | Fields or sections | Meaning |
|---|---|---|
| `report_facts.json` | Per-stage counts, result categories, top candidates, and thresholds | Machine-readable source of the report's numbers. |
| `report.md` / `report.html` | Main findings, biological interpretation, limitations, and next steps | Human-readable interpretation that separates observations from hypotheses. |
| `agent_metadata.json` | Agent status and model metadata | Records whether an optional plain-language interpretation was added. |

Any optional language-model interpretation is based only on the compact facts file; raw genotypes and full association tables remain local.

### Implementation notes

Stage 9 is an optional, separate pipeline command. Install its dependencies only when
reporting is wanted:

```bash
pip install -r claude-first-prototype/requirements-stage9.txt

python claude-first-prototype/pipeline/stage9_report.py \
  --config claude-first-prototype/stage5_output/config.yaml \
  --config claude-first-prototype/stage8_output/config.yaml \
  --out-dir claude-first-prototype/stage9_output \
  --agent required
```

Set `OPENAI_API_KEY` before using `--agent required`. `--agent auto` uses the agent
when a key is available, while `--agent off` produces only the deterministic report.
Stages 1-8 do not import the OpenAI SDK, install Stage 9 dependencies, or invoke Stage 9.
The optional agent interpretation should always be read as a summary of pipeline
results, not as independent biological evidence.

## Software

The implemented stages are Python-based. Stage 5 uses SciPy. Stage 0 is expected to use `bcftools` and `truvari` when implemented.
