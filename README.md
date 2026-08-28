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
gendiff old.bam new.bam --reference reference.fa --igv-batch differences.igv.batch
gendiff old.bam new.bam --threads 8 --progress
gendiff old.vcf.gz new.bcf --json
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
record matching to report identical, modified, added, and removed records; use
`--temp-dir` to select scratch storage. Indexed alignment files can be divided by
contig when more than two threads are requested. Other inputs still scan both
files concurrently.

Reference-aware VCF normalization follows bcftools semantics and requires an
indexed FASTA. HTML reports are self-contained, while SVG summaries can be used
directly in documents and presentations. IGV output is a batch file that loads
both inputs and visits example changed loci.

Some Linux distributions provide an unrelated `/usr/bin/gendiff`; install GenDiff
in a virtual environment or invoke it as `python -m gendiff`. GenDiff is not a
substitute for assay validation.

## License

MIT
