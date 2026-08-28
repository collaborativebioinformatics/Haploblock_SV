# Stage 0: Merge and collapse a cohort of SV VCFs

[`pipeline/stage_0_mergingvcf.py`](pipeline/stage_0_mergingvcf.py) accepts a
manifest of single-sample VCF/VCF.GZ files, prepares them for cohort merging,
and writes a multi-sample VCF plus `truvari`-collapsed outputs.

## Requirements

- Python 3.9 or newer
- `bcftools`
- `truvari`
- An indexed reference FASTA (`.fai`) when `--reference` is supplied

Each manifest line is one local VCF path. Blank lines and `#` comments are
ignored. At least two files are required, and each input must contain exactly
one sample. Inputs must use the same contig names and lengths unless
`--force-mismatched-references` is explicitly supplied.

## Example

```bash
python3 pipeline/stage_0_mergingvcf.py \
  --vcf-list giab_data/list_vcf.txt \
  --out-dir stage0_vcf_merge_output \
  --reference GRCh38.fa \
  --threads 4
```

Use `--allow-duplicate-sample-names` only when duplicate sample names are
intentional; it passes `--force-samples` to `bcftools merge`. Tune matching
with repeated `--truvari-arg` options, for example
`--truvari-arg=--refdist --truvari-arg=500`.

## Processing

1. Validate the manifest, sample names, and reference contigs.
2. Run `bcftools sort` on every input.
3. Assign cohort-unique IDs with `bcftools annotate`, then index each staged
   input.
4. Run `bcftools merge -m id` so symbolic SVs are not incorrectly joined by
   position alone.
5. Run `truvari collapse`.
6. Sort **both** `truvari` outputs before creating tabix indexes. This is
   required because `truvari` does not guarantee coordinate order.

## Outputs

```text
stage0_vcf_merge_output/
├── staged/
│   ├── sample0000.sorted.vcf.gz
│   ├── sample0000.staged.vcf.gz
│   ├── sample0000.staged.vcf.gz.tbi
│   └── ...
├── merged.exact.vcf.gz
├── merged.exact.vcf.gz.tbi
├── merged.collapsed.vcf.gz
├── merged.collapsed.vcf.gz.tbi
├── merged.collapsed_removed.vcf.gz
├── merged.collapsed_removed.vcf.gz.tbi
└── run.json
```

`merged.exact.vcf.gz` preserves every post-merge record. The collapsed VCF
contains records retained by `truvari`; `merged.collapsed_removed.vcf.gz`
contains records collapsed into another record. The final VCFs are sorted and
tabix-indexed for direct use by Stage 1.

`run.json` records input paths, options, staged files, and output paths. A
non-zero exit status means the run is incomplete and its outputs should not be
used without inspection.

## Scope

Stage 0 performs representation-level merging and collapsing only. It does
not run kanpig, apply analysis-specific size filters, or silently remove
`IMPRECISE` records. Downstream stages decide which records to analyze.
