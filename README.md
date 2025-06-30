## Usage

```sh
Checks integrity and consistency of FASTQ files.

Use --paired or --single for each sample. These flags can be used multiple times.

By default, the tool will exit immediately after the first error is found. Use --continue-on-error to check all files regardless of errors.

Usage: grz-check fastq [OPTIONS] --output <OUTPUT>

Options:
      --paired <FQ1_PATH> <FQ2_PATH> <FQ1_LEN> <FQ2_LEN>
          A paired-end sample. Provide FQ1, FQ2, FQ1 length, and FQ2 length. Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check. Example: --paired fq1.fastq.gz fq2.fastq.gz 150 151

      --single <FQ_PATH> <READ_LEN>
          A single-end sample. Provide the file path and read length. Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check. Example: --single sample.fastq.gz 0

  -o, --output <OUTPUT>
          Path to write the output TSV report

      --continue-on-error
          Continue processing all files even if an error is found

      --threads <THREADS>
          Number of threads to use for processing

  -h, --help
          Print help (see a summary with '-h')
```

## Example

```sh
# --paired R1 R2 read_length_R1 read_length_R2
# --single R1 read_length_R1
grz-check --show-progress fastq --output report.tsv --paired path/to/sample__R1.fastq.gz path/to/sample_R2.fastq.gz 150 150 --single path/to/sample.fastq.gz 151
```
