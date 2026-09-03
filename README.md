# SemantiSeq

Semantic comparison and regression testing for genomics files.

SemantiSeq determines whether two files contain the same logical records, even when
record order, compression, or provenance metadata differ. It supports BAM/CRAM
and VCF/BCF.

```text
$ semantiseq sample-a.merged.bam sample-b.merged.bam --explain
Equivalent: no
Type: alignment
Relationship: different datasets
Identity overlap: 0.0%
Content similarity: 0.0%
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

![Example SemantiSeq comparison](assets/example-output.svg)

## Install

Python 3.9 or newer is required.

```bash
pip install semantiseq
```

For local development:

```bash
git clone https://github.com/AlexanderM-M/semantiseq.git
cd semantiseq
python -m pip install -e '.[test]'
pytest
```

## Usage

```bash
semantiseq OLD NEW
semantiseq old.bam new.bam --explain
semantiseq old.bam new.bam --profile mapping --ignore-tag MD
semantiseq old.cram new.cram --reference reference.fa
semantiseq old.vcf.gz new.bcf --normalize --reference reference.fa
semantiseq old.bam new.bam --html comparison.html
semantiseq old.bam new.bam --svg comparison.svg
semantiseq old.bam new.bam --write-diff diff-records
semantiseq old.bam new.bam --tracks genomic-tracks
semantiseq old.bam new.bam --difference-table differences.tsv.gz
semantiseq old.bam new.bam --region chr1:100000-200000
semantiseq old.bam new.bam --regions targets.bed --exclude-regions blacklist.bed
semantiseq old.bam new.bam --reference reference.fa --igv-batch differences.igv.batch
semantiseq old.bam new.bam --threads 8 --progress
semantiseq old.bam new.bam --cache /scratch/semantiseq-cache
semantiseq old.bam new.bam --reproducibility comparison.json
semantiseq old.vcf.gz new.bcf --json
semantiseq batch comparisons.tsv --html regression.html --junit regression.xml
semantiseq manifest baseline-run candidate-run --output comparisons.tsv
semantiseq matrix run/*.bam --html cohort.html --cache /scratch/semantiseq-cache
semantiseq identity old.vcf.gz new.vcf.gz
semantiseq identity old.bam new.bam --sites fingerprint-sites.vcf.gz
```

Exit status is `0` when inputs are semantically equivalent, `1` when they differ,
and `2` when comparison cannot be completed. JSON output is intended for CI and
workflow integration.

SemantiSeq compares all logical record fields and the structural parts of each
header. Record ordering, file encoding, and non-structural header metadata do not
affect equivalence. The `strict`, `core`, `mapping`, `calls`, and `genotypes`
profiles provide progressively focused comparisons for their supported formats.
Input names are inferred from filenames and can be overridden with `--name-a` and
`--name-b`.

The default fingerprint pass uses bounded memory. `--explain` performs disk-backed
record matching and reports both changed fields and transitions such as MAPQ,
mapping status, flags, filters, and genotypes. Use `--temp-dir` to select scratch
storage. With indexed alignment files, both scanning and detailed matching are
divided by contig when more than two threads are requested.

`--cache DIR` reuses validated semantic scans and disk-backed matching data. The
key includes the input path, size, nanosecond modification time, SemantiSeq version,
reference metadata, profile, ignored fields, and region filters. Completed input
scans survive an interrupted comparison and are resumed on the next invocation.
Caches created for detailed comparisons contain record summaries and should be
protected like the source data.

Limit a comparison with repeatable `--region` values or a BED file passed to
`--regions`; `--exclude-regions` removes unwanted intervals. Coordinates passed
to `--region` are one-based and inclusive, while BED coordinates follow the BED
standard. `--difference-table` writes every modified or sample-only record as a
streamed TSV ledger; use a `.gz` suffix for compression.

`--write-diff DIR` writes four ordinary BAM or VCF files: records found only in
each input and the before/after versions of modified records. A small manifest
records the labels, comparison profile, filenames, and record counts. These files
can be passed directly to samtools, bcftools, IGV, or downstream workflows.

`--tracks DIR` writes affected loci, raw and normalized difference bedGraphs, a
log2 coverage-ratio bedGraph, and a per-contig TSV. These files work directly
with IGV, bedtools, and other genome browsers.

Reference-aware VCF normalization follows bcftools semantics and requires an
indexed FASTA. HTML reports are self-contained and lead with three concise
findings, a proportional record overview, significant distribution shifts, and
a normalized genome-wide difference and coverage map. A MAPQ or genotype
transition matrix appears when applicable; detailed tables stay collapsed. SVG
summaries can be used directly in documents and presentations. IGV output is a
batch file that loads both inputs and visits example changed loci.

Header diagnostics identify reference, contig, read-group, sample, INFO, FORMAT,
and FILTER definition changes. Record differences are attributed to BAM read
groups or VCF samples when possible. `--reproducibility FILE` writes input and
reference SHA-256 checksums, the exact command, settings, result, version, and
runtime.

## Cohorts and sample identity

`semantiseq matrix` scans each file once and creates a clustered semantic-similarity
matrix, with likely outliers highlighted. It accepts the normal profile,
normalization, reference, and ignored-field options.

`semantiseq identity` compares genotype fingerprints and reports the best sample
matches. VCF/BCF inputs work directly. Indexed BAM/CRAM inputs require a
biallelic SNP panel through `--sites`; CRAM may additionally require
`--reference`. Low-call results are reported as insufficient rather than treated
as a match. Identity checking is a swap guard, not a contamination estimate.

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
semantiseq batch comparisons.tsv \
  --threads 16 \
  --html regression.html \
  --json regression.json \
  --junit regression.xml \
  --multiqc semantiseq_mqc.json \
  --write-diff differences
```

By default, modified records, records found in only one input, and structural
header differences fail a comparison. Controlled changes can be accepted with
`--max-modified`, `--max-only`, `--min-overlap`, and
`--allow-structural-differences`. `--fail-transition genotype`, for example,
continues to reject genotype transitions even when some modified records are
allowed. Optional manifest columns are `profile`, `reference`, and `normalize`;
relative paths are resolved from the manifest directory.

For stage-specific acceptance rules, pass a versioned JSON policy:

```bash
semantiseq batch comparisons.tsv --policy examples/policy.json --json semantiseq-batch.json
```

Rules can match sample, stage, file kind, and profile. They support absolute and
relative limits, allowed or forbidden changed fields and transitions, and
error, warning, or informational severity. Reports record the matched rule and
the decision trace. Unapproved transitions fail by default when a policy file is
used, and later matching rules override earlier settings.

To create the manifest from matching directory trees, use `semantiseq manifest`.
It pairs relative paths exactly and stops if either tree has missing or ambiguous
files. Add `--dry-run` to inspect the TSV without writing it.

Compact baselines preserve semantic hashes and counts without storing read names,
loci, sequences, or whole input files:

```bash
semantiseq baseline create comparisons.tsv --output pipeline-baseline.semantiseq.json
semantiseq baseline check pipeline-baseline.semantiseq.json comparisons.tsv
```

Baseline checks detect changes to logical fields and structural headers. Keep the
original baseline files when record-level diff files or transition explanations
may be needed later.

The aggregate JSON is intended for automation, JUnit integrates with CI systems,
and filenames ending in `_mqc.json` are discovered by MultiQC as custom content.
Installing `semantiseq[multiqc]` also provides a native MultiQC module; name the
aggregate JSON `semantiseq-batch.json` so it is discovered automatically.

SemantiSeq can also be invoked as `python -m semantiseq`. It is not a substitute
for assay validation.

## License

MIT
