# Methods documentation

This directory contains detailed methods documentation for the prototype pipeline. The user-facing [README](../README.md) explains the biological questions and how to interpret the main outputs; these documents describe inputs, data contracts, calculations, parameters, and limitations.

| Stage | Methods document | Required upstream output |
|---|---|---|
| 0 | [Cohort SV merging](stage0_cohort_sv_merging.md) | Single-sample SV VCFs (planned) |
| 1 | [Cluster-aware preprocessing](stage1_cluster_aware_preprocessing.md) | Cohort SV VCF |
| 2 | [Boundary classification](stage2_boundary_classification.md) | Stage 1 configuration |
| 4 | [Population frequency classification](stage4_population_frequency.md) | Stage 1 configuration |
| 5 | [SV-type enrichment](stage5_sv_type_enrichment.md) | Stage 1 configuration |
| 6 | [Population-conditioned association](stage6_cluster_association.md) | Stage 4 configuration |
| 7 | [Hash representation audit](stage7_hash_representation.md) | Stage 6 configuration |
| 8 | [Candidate annotation](stage8_candidate_annotation.md) | Stage 7 configuration |
| 9 | [Reporting](stage9_reporting.md) | Stage 5 and Stage 8 configurations |

Stages are numbered to match their scripts. There is no Stage 3 analysis in the prototype.
