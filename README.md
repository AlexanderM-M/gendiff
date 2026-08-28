# GenDiff

Semantic comparison for genomics files.

GenDiff determines whether two files contain the same logical records, even when
record order, compression, or provenance metadata differ. It supports BAM/CRAM
and VCF/BCF.

```console
$ gendiff old.bam new.bam
Equivalent: yes
Type: alignment
Records: 84122931 / 84122931
Logical records: same
Structural header: same
Metadata header: different (informational)
```

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
gendiff old.cram new.cram --reference reference.fa
gendiff old.vcf.gz new.bcf --json
gendiff old.bam new.bam --threads 8
```

Exit status is `0` when inputs are semantically equivalent, `1` when they differ,
and `2` when comparison cannot be completed. JSON output is intended for CI and
workflow integration.

GenDiff compares all logical record fields and the structural parts of each
header. Record ordering, file encoding, and non-structural header metadata do not
affect equivalence. It reports which field groups differ without loading complete
files into memory. Two worker threads are used by default; increase `--threads`
when comparing large files on a machine or allocated cluster job with more cores.

GenDiff does not currently normalize alternative variant representations. Run
variant normalization before comparing files when equivalent variants may be
encoded differently. It is not a substitute for assay validation.

## License

MIT
