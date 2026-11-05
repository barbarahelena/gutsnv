# gutsnv: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - [unreleased]

Initial version, written with Nextflow static typing (≥ 26.10.5).

### `Added`

- Read preprocessing: run merging, FastQC, fastp, Bowtie2 host removal, seqkit read counts (each optional)
- GT-Pro metagenotyping with decoding into per-site allele counts and species bins
- Per-species allele-frequency matrices, strain quasi-phaseability and per-cohort SNV/species QC
- Per-species strain distance, kinship matrix and PCoA
- Species abundance: GT-Pro-based estimate (always) and CoverM (optional, `--run_coverm`)
- Summary tables and MultiQC report
- Typed modules: fastqc, fastp, cat/fastq, bowtie2/build, bowtie2/align (as host removal), seqkit/stats, coverm/genome and multiqc
- GT-Pro container on Docker Hub (`barbarahelena/gt-pro:1.0.1`)
- nf-test stub tests (full run and all optional steps skipped)
