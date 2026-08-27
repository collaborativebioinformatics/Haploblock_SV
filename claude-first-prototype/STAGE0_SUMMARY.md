# Stage 0: Merging and Collapsing VCF Files

[`pipeline/stage_0_mergingvcf.py`](../pipeline/stage_0_mergingvcf.py) stages
single-sample structural-variant VCFs, merges them with `bcftools`, and
collapses nearby equivalent calls with `truvari`.

## Requirements

The following executables must be available on `PATH`:

- `bcftools`
- `truvari`

The file supplied to `--vcf-list` must contain one VCF or VCF.GZ path per
line. Blank lines and lines beginning with `#` are ignored.

## Example command

```bash
python3 stage_0_mergingvcf.py \
  --vcf-list giab_data/list_vcf.txt \
  --out-dir stage0_vcf_merge_output \
  --allow-duplicate-sample-names
```

`--allow-duplicate-sample-names` passes `--force-samples` to `bcftools merge`.
Without this option, duplicate sample names are rejected before merging.

## Workflow and outputs

For each input VCF, the script:

1. Sorts the file with `bcftools sort`.
2. Assigns a unique record ID with `bcftools annotate`.
3. Indexes the staged VCF.

It then:

1. Merges staged files with `bcftools merge -m id`.
2. Indexes the exact merge.
3. Runs `truvari collapse`.
4. Indexes the collapsed and removed-record VCFs.

Expected output files are:

```text
/beegfs/l233n428/hack/stage0_vcf_merge_output3/
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

The exact merge preserves records before approximate matching. The collapsed
VCF contains the records retained by `truvari`; the
`merged.collapsed_removed.vcf.gz` file contains records collapsed into another
record.

## Interpreting the reported run

The log shows that sorting, annotation, staging indexes, the exact merge, its
index, and `truvari collapse` completed successfully. `truvari` also reported:

- 125,210 input variants
- 125,163 variants written to the collapsed output
- 47 variants collapsed into 46 variants

However, the run did **not** complete successfully. The final command failed:

```text
[E::hts_idx_push] Chromosome blocks not continuous
index: failed to create an index for ".../merged.collapsed.vcf.gz"
ERROR: Command failed with exit code 255: bcftools
```

This means `merged.collapsed.vcf.gz` was written, but its tabix index was not.
The output should not be treated as a complete Stage 0 result because
downstream tools commonly require the `.tbi` index. The script should also
have exited with status `1`; check this with:

```bash
echo $?
```

The error indicates that the collapsed VCF is not in the coordinate order
required by tabix. `truvari collapse` can change the order of retained
records, so the collapsed output should be sorted before indexing. A
workaround is:

```bash
out=/beegfs/l233n428/hack/stage0_vcf_merge_output3

mv "$out/merged.collapsed.vcf.gz" "$out/merged.collapsed.unsorted.vcf.gz"
bcftools sort -Oz \
  -o "$out/merged.collapsed.vcf.gz" \
  "$out/merged.collapsed.unsorted.vcf.gz"
bcftools index --tbi "$out/merged.collapsed.vcf.gz"
```

The removed-record output should be checked and sorted similarly if it also
fails indexing:

```bash
bcftools sort -Oz \
  -o "$out/merged.collapsed_removed.sorted.vcf.gz" \
  "$out/merged.collapsed_removed.vcf.gz"
bcftools index --tbi "$out/merged.collapsed_removed.sorted.vcf.gz"
```

After sorting, verify that the index exists and that the VCF can be queried:

```bash
test -f "$out/merged.collapsed.vcf.gz.tbi"
bcftools index -n "$out/merged.collapsed.vcf.gz"
```

The durable fix is to update `stage_0_mergingvcf.py` so it runs
`bcftools sort` on each `truvari` output before calling `bcftools index`.
