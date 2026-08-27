> Structural Variants Hackathon at Baylor College of Medicine, August 25-28, 2026

# Haploblock_SV

# Overview

We aim to reveal how structural variation is organized within and around haploblocks, and what that says about population structure. This is a Follow up on previous projects: [haploblocks.org](https://haploblocks.org) and [data.haploblocks.org](https://data.haploblocks.org)

Specifically, the workflow will: 
- Identify structural variants (SVs) including deletions, duplications, inversions and insertions that occur within haploblocks, and determine whether any SV-breakpoints overlap or occur near haploblock boundaries.
- Classifies SVs as common or population-specific using population allele frequency data, creating seperate datasets for downstream comparitive analysis. 
- Assess haploblock enrichment for specifc SV types by testing whether certain haploblocks contain significanlty higher or lower numbers of particular SV classess than expected based on haploblock size and overall SV distribution
- Evaluate relationship between SV patterns and population structure by comparing SV classifiers and distributions with population cluster derived from 1000 Genome small-variant hashes [data.haploblocks.org] 
- Develp a reusable and scalable analysis pipeline that automates the SV annotation, classification, enrichment testing, and population strucutre comparisons, allowing the workflow applied efficiently to additional datasets and future studies. 


## Contributors
- Jędrzej Kubica jedrzej.kubica@univ-grenoble-alpes.fr 
- Lynn Ly lynn.ly@nanoporetech.com
- Maria Fernanda Cardenas maria.cardenas@stjude.org
- Linh Nguyen nguyen.linh.1010@ku.edu


## Flowchart

![flowchart_haploblock_SV_upd.png](flowchart_haploblock_SV_upd.png)

## Data

- genomic hashes and clusters: [data.haploblocks.org](https://data.haploblocks.org)
- 1000 Genomes HGSVC: [https://www.internationalgenome.org/human-genome-structural-variation-consortium/](https://www.internationalgenome.org/human-genome-structural-variation-consortium)
- Long read 1KGP SV calls: [https://s3.amazonaws.com/1000g-ont/index.html?prefix=PROCESSED_DATA/ALIGNED_TO_HG38/SNIFFLES_v2.6.2/](https://s3.amazonaws.com/1000g-ont/index.html?prefix=PROCESSED_DATA/ALIGNED_TO_HG38/SNIFFLES_v2.6.2/)
- 1KGP haplotype-resolved SVs [https://www.nature.com/articles/s41467-018-08148-z](https://www.nature.com/articles/s41467-018-08148-z): Study ID nstd152 (Chaisson et al. 2019)

## Pipeline Overview 

| Stage | Name | Purpose |
|---|---|---|
| Stage 0 | Data ingestion and harmonization | Download and reconcile VCF files |
| Stage 1 | Data preparation | Prepare and standardize input variant and haplotype data |
| Stage 2 | SV annotation | Annotate SVs relative to haploblocks |
| Stage 3 | Boundary enrichment | Removed |
| Stage 4 | Haploblock classification | Define common vs population-specific haploblocks |
| Stage 5 | SV type enrichment | Test enrichment of SV types across haploblock classes |
| Stage 6 | Cluster comparison | Compare haploblock/SV patterns with previously defined population clusters |


## Development status 

The pipeline is currently under active development 

- [x] Repository setup
- [x] Define pipeline architecture
- [x] Stage 0 — Data ingestion and harmonization
- [x] Stage 1 — Data preparation
- [x] Stage 2 — SV annotation
- [ ] Stage 3 — Boundary enrichment
- [ ] Stage 4 — Haploblock classification
- [ ] Stage 5 — SV type enrichment
- [ ] Stage 6 — Population cluster comparison


# Results

# Future steps

- Generate Region-level FASTA including SNP and SV information to regenerate Clustering Haploblock with MMseqs2 and to confirm that the genomic hashes don't have any deviation with the ones originated previously 


# References
