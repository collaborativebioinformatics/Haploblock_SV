> Structural Variants Hackathon at Baylor College of Medicine, August 25-28, 2026

# Haploblock_SV

# Quickstart
Install dependencies:
```bash
pip install -r claude-first-prototype/requirements.txt
```

Run the pipeline:
```bash
python claude-first-prototype/pipeline/stage1_cluster_aware.py
```

The current workflow starts from a previously merged VCF and writes chromosome-specific SV-to-cluster association tables into `stage1_output/`.
See [claude-first-prototype/README.md](claude-first-prototype/README.md) for output fields, interpretation, QC, and the optional boundary-position analysis.
Stage 0 is currently a placeholder to download single-sample VCFs from 1kgp and merge them into a cohort VCF.

# Overview

We aim to reveal how structural variation is organized within and around haploblocks, and what that may tell us about haplotype and population structure.

This project follows [haploblocks.org](https://haploblocks.org) and [data.haploblocks.org](https://data.haploblocks.org), which determined haploblock boundaries and clustered haplotypes within each block using 1000 Genomes small variants. Those results did not incorporate structural variants, so they cannot directly identify which existing haplotype clusters carry a particular SV.

So far implemented: cluster-aware preprocessing that combines cohort SV genotypes with those existing cluster memberships and probabilistically assigns SVs to one or more clusters. This is necessary because the VCF and haploblock clusters were phased independently, so their `hap0` and `hap1` labels cannot be matched directly.

The wider pipeline is exploratory rather than a fixed production analysis. It is intended to test biological expectations about SVs and haploblocks—for example, whether SVs consistently mark particular haplotype clusters, whether some blocks or clusters contain unusual SV patterns, and whether those patterns reflect population structure. The downstream analyses will evolve as these assumptions are evaluated.

- Infer which haplotype clusters within each haploblock carry each SV (currently DEL and INS)
- Separately, if useful, check whether SVs fall near or across haploblock boundaries
- Assess Haploblock enrichment for specific SV types by testing whether certain haploblocks contain significantly higher or lower number of particular SV classes.
- Classify structural variants (SVs) as common or population-specific based on population frequency data
- Check if these SVs correlate with 1000 Genomes population clusters
- Wrap it up in a re-usable and scalable pipeline


## Contributors
- Jędrzej Kubica jedrzej.kubica@univ-grenoble-alpes.fr 
- Lynn Ly lynn.ly@nanoporetech.com
- Maria Fernanda Cardenas maria.cardenas@stjude.org
- Linh Nguyen nguyen.linh.1010@ku.edu

# Workflow

![Cluster-aware Haploblock SV workflow](flowchart_haploblock_SV_upd.png)

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
- [x] Stage 4 — Common vs population-specific SV classification
- [ ] Stage 5 — SV type enrichment
- [ ] Stage 6 — Population cluster comparison


# Results

# Future steps

- Create region-level FASTA files containing SNPs and SVs, reclusters the haploblocks using MMseqs2, and verify that the resulting clusters match the original ones. 


# References
https://osf.io/preprints/biohackrxiv/xhkc3_v1
