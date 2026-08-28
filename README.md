# GenDiff

Semantic comparison for genomics files.

GenDiff determines whether two files contain the same logical records, even when
record order, compression, or provenance metadata differ. It supports BAM/CRAM
and VCF/BCF.

```text
$ gendiff sample-a.merged.bam sample-b.merged.bam --explain
Equivalent: no
Type: alignment
Relationship: different datasets
Identity overlap: 0.0%
Inputs:
  Sample A: sample-a.merged.bam (9,398 records)
  Sample B: sample-b.merged.bam (147,213 records)
Logical records: different
Structural header: different
Metadata header: different (informational)
Record differences:
  Identical: 0
  Modified: 0
  Only in Sample A: 9,398
  Only in Sample B: 147,213
```

![Example GenDiff comparison](assets/example-output.svg)

## Install

Python 3.9 or newer is required.

```bash
python -m pip install git+https://github.com/AlexanderM-M/gendiff.git
```

For local development:

```bash
git clone https://github.com/AlexanderM-M/gendiff.git
cd gendiff
python -m pip install -e '.[test]'
pytest
```

## Usage

```bash
gendiff OLD NEW
gendiff old.bam new.bam --explain
gendiff old.bam new.bam --profile mapping --ignore-tag MD
gendiff old.cram new.cram --reference reference.fa
gendiff old.vcf.gz new.bcf --normalize --reference reference.fa
gendiff old.bam new.bam --html comparison.html
gendiff old.bam new.bam --svg comparison.svg
gendiff old.bam new.bam --write-diff diff-records
gendiff old.bam new.bam --reference reference.fa --igv-batch differences.igv.batch
gendiff old.bam new.bam --threads 8 --progress
gendiff old.vcf.gz new.bcf --json
gendiff batch comparisons.tsv --html regression.html --junit regression.xml
```

Exit status is `0` when inputs are semantically equivalent, `1` when they differ,
and `2` when comparison cannot be completed. JSON output is intended for CI and
workflow integration.

GenDiff compares all logical record fields and the structural parts of each
header. Record ordering, file encoding, and non-structural header metadata do not
affect equivalence. The `strict`, `core`, `mapping`, `calls`, and `genotypes`
profiles provide progressively focused comparisons for their supported formats.
Input names are inferred from filenames and can be overridden with `--name-a` and
`--name-b`.

The default fingerprint pass uses bounded memory. `--explain` performs disk-backed
record matching and reports both changed fields and transitions such as MAPQ,
mapping status, flags, filters, and genotypes. Use `--temp-dir` to select scratch
storage. Indexed alignment files can be divided by contig when more than two
threads are requested. Other inputs still scan both files concurrently.

`--write-diff DIR` writes four ordinary BAM or VCF files: records found only in
each input and the before/after versions of modified records. A small manifest
records the labels, comparison profile, filenames, and record counts. These files
can be passed directly to samtools, bcftools, IGV, or downstream workflows.

Reference-aware VCF normalization follows bcftools semantics and requires an
indexed FASTA. HTML reports are self-contained, while SVG summaries can be used
directly in documents and presentations. IGV output is a batch file that loads
both inputs and visits example changed loci.

## Pipeline regression

Batch mode compares multiple baseline and candidate outputs. The manifest is a
tab-separated file; row order defines pipeline-stage order for earliest-divergence
reporting.

```text
sample  stage       before                 after
Sample A alignment  baseline/A.bam         candidate/A.bam
Sample A variants   baseline/A.vcf.gz       candidate/A.vcf.gz
Sample B variants   baseline/B.vcf.gz       candidate/B.vcf.gz
```

```bash
gendiff batch comparisons.tsv \
  --threads 16 \
  --html regression.html \
  --json regression.json \
  --junit regression.xml \
  --multiqc gendiff_mqc.json \
  --write-diff differences
```

By default, modified records, records found in only one input, and structural
header differences fail a comparison. Controlled changes can be accepted with
`--max-modified`, `--max-only`, `--min-overlap`, and
`--allow-structural-differences`. `--fail-transition genotype`, for example,
continues to reject genotype transitions even when some modified records are
allowed. Optional manifest columns are `profile`, `reference`, and `normalize`;
relative paths are resolved from the manifest directory.

The aggregate JSON is intended for automation, JUnit integrates with CI systems,
and filenames ending in `_mqc.json` are discovered by MultiQC as custom content.

Some Linux distributions provide an unrelated `/usr/bin/gendiff`; install GenDiff
in a virtual environment or invoke it as `python -m gendiff`. GenDiff is not a
substitute for assay validation.

## License

MIT
