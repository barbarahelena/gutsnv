# gutsnv

[![tests](https://github.com/barbarahelena/gutsnv/actions/workflows/nf-test.yml/badge.svg)](https://github.com/barbarahelena/gutsnv/actions/workflows/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A526.10.5-0DC09D?logo=nextflow&logoColor=white)](https://www.nextflow.io/)

Metagenomic SNV genotyping of gut microbes with [GT-Pro](https://github.com/zjshi/gt-pro), strain quasi-phaseability, per-cohort SNV QC and per-species strain kinship matrices: the inputs for a microbial SNV association study (SNV-MWAS).

The workflow follows the structure of Wang et al., _Nature_ 2026, "Microbiome-wide PheWAS links gut microbial SNVs to human health and exposures" ([code](https://github.com/GRONINGEN-MICROBIOME-CENTRE/SNV-MWAS)), except for the preprocessing steps of the sequencing data (fastp and Bowtie2 instead of Kneaddata). Since the downstream association analysis is study-specific, these analyses are not part of this workflow.

## Pipeline overview

1. **Read preprocessing** (each step optional)
   - Merge runs/lanes per sample
   - Raw read QC ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
   - Adapter and quality trimming ([fastp](https://github.com/OpenGene/fastp)), with settings equivalent to the paper's Trimmomatic settings
   - Host-read removal ([Bowtie2](https://github.com/BenLangmead/bowtie2))
   - Read counts of the reads entering genotyping ([seqkit](https://github.com/shenwei356/seqkit)), the sequencing-depth covariate
2. **Metagenotyping**: GT-Pro detects known SNVs of gut species by exact k-mer matching. Raw output is decoded into reference/alternative allele counts per site.
3. **Per-species genotypes**, computed in parallel species bins:
   - alternative-allele frequency per sample and SNV
   - strain **quasi-phaseability** per sample and species: whether a sample is dominated by a single strain
   - per-cohort SNV QC (missingness, MAF) and species QC (minimum numbers of SNVs and samples)
4. **Strain kinship** per species: scaled Manhattan genetic distance, kinship matrix (1 − distance) and PCoA
5. **Species abundance**: a GT-Pro-based abundance estimate per species (always), and optionally [CoverM](https://github.com/wwood/CoverM) read mapping to one reference genome per species (`--run_coverm`)
6. **Summary** tables of the samples, quasi-phaseable samples and SNVs available per species and cohort, plus a [MultiQC](https://multiqc.info) report

## Requirements

- Nextflow **≥ 26.10.5**. The pipeline uses static typing: typed params, records, typed processes and workflow outputs.
- A container engine (Docker, Singularity/Apptainer) or an institutional config profile (e.g. `-profile snellius`).
- The GT-Pro database: the `<name>_optimized_db_kmer_index.bin` and `<name>_optimized_db_snps.bin` files plus the SNP dictionary (`variants_main.covered.hq.snp_dict.tsv`). See the [GT-Pro README](https://github.com/zjshi/gt-pro). The public database is the 881-species build (`20190723_881species`).
- For host removal: a host genome FASTA or a Bowtie2 index.
- For abundance: a directory with one reference genome FASTA per species.

## Usage

Prepare a samplesheet with one row per sequencing run:

```csv
sample,fastq_1,fastq_2,cohort,subject,timepoint
SAMPLE_A,/path/SAMPLE_A_L001_R1.fastq.gz,/path/SAMPLE_A_L001_R2.fastq.gz,cohort1,subject1,t1
SAMPLE_A,/path/SAMPLE_A_L002_R1.fastq.gz,/path/SAMPLE_A_L002_R2.fastq.gz,cohort1,subject1,t1
SAMPLE_B,/path/SAMPLE_B.fastq.gz,,cohort2,subject2,
```

`cohort` defines the groups in which SNV QC is applied (per cohort, as in the paper). `subject` and `timepoint` are carried through for repeated-measures designs.

```bash
nextflow run barbarahelena/gutsnv \
    -profile <docker/singularity/institute> \
    --input samplesheet.csv \
    --gtpro_db /path/to/gtpro_db \
    --gtpro_snp_dict /path/to/gtpro_db/variants_main.covered.hq.snp_dict.tsv \
    --host_index /path/to/bowtie2/GRCh38 \
    --outdir results
```

For reads that are already trimmed and host-depleted, add `--skip_qc --skip_host_removal`. See [docs/usage.md](docs/usage.md) for all parameters and [docs/output.md](docs/output.md) for the outputs.

## Credits

The analysis is based on the SNV-MWAS pipeline of Wang et al. (2026), _Microbiome-wide PheWAS links gut microbial SNVs to human health and exposures_, Nature, doi: [10.1038/s41586-026-11220-0](https://doi.org/10.1038/s41586-026-11220-0) ([GitHub](https://github.com/GRONINGEN-MICROBIOME-CENTRE/SNV-MWAS)).

The pipeline was built using the [nf-core](https://nf-co.re) template (Ewels et al. 2020, Nat Biotechnol, doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)).

## Citations

A list of the tools used by the pipeline is in [CITATIONS.md](CITATIONS.md).
