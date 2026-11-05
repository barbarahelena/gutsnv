# gutsnv: Usage

## Samplesheet

`--input` is a CSV file with a header row. Columns are matched by name:

| Column      | Required | Description                                                                          |
| ----------- | -------- | ------------------------------------------------------------------------------------ |
| `sample`    | yes      | Sample ID. Rows with the same ID are runs/lanes of one sample and are concatenated.  |
| `fastq_1`   | yes      | Gzipped FASTQ (R1).                                                                  |
| `fastq_2`   | no       | Gzipped FASTQ (R2). Leave empty for single-end data.                                 |
| `cohort`    | no       | Group in which SNV QC is applied. Default: `all`.                                    |
| `subject`   | no       | Person ID, e.g. for repeated measures. Default: the sample ID.                       |
| `timepoint` | no       | Timepoint label.                                                                     |

The samplesheet is read and validated by Nextflow itself against the `Sample` record type (`types/records.nf`).

## Reference data

| Parameter               | Description                                                                                                               |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `--gtpro_db`            | Directory with the optimised GT-Pro database files `<gtpro_db_name>_optimized_db_kmer_index.bin` and `..._snps.bin`.      |
| `--gtpro_db_name`       | Database prefix. Default: `20190723_881species`.                                                                         |
| `--gtpro_snp_dict`      | GT-Pro SNP dictionary (`variants_main.covered.hq.snp_dict.tsv`), plain or gzipped.                                        |
| `--host_index`          | Directory with a Bowtie2 index of the host genome; or                                                                     |
| `--host_fasta`          | Host genome FASTA (an index is built).                                                                                    |
| `--abundance_genomes`   | Directory with one reference genome FASTA per species; required with `--run_coverm`.                                                       |
| `--genome2species`      | Optional TSV without header, `genome<TAB>species_id`, to name abundance columns by GT-Pro species ID.                     |

## Parameters and defaults

All parameters are declared, with types and defaults, in the `params` block of `main.nf`. `nextflow_schema.json` is generated from that block (`python .github/scripts/params_to_schema.py main.nf nextflow_schema.json`).

### Read preprocessing

| Parameter             | Default | Description                                     |
| --------------------- | ------- | ----------------------------------------------- |
| `--skip_fastqc`       | false   | Skip FastQC on raw reads.                       |
| `--skip_qc`           | false   | Skip adapter/quality trimming with fastp.       |
| `--skip_host_removal` | false   | Skip host-read removal.                         |
| `--save_clean_reads`  | false   | Publish the reads that enter genotyping.        |

Trimming uses `--cut_right --cut_window_size 4 --cut_mean_quality 20 --length_required 50`, equivalent to Trimmomatic `SLIDINGWINDOW:4:20 MINLEN:50`. Host removal uses Bowtie2 `--very-sensitive-local` (the KneadData default). Both can be changed via `ext.args` in a custom config.

### Species abundance

| Parameter                       | Default | Description                                                                             |
| ------------------------------- | ------- | --------------------------------------------------------------------------------------- |
| `--run_coverm`                  | false   | Also estimate relative abundance with CoverM (read mapping; computationally expensive). |
| `--abundance_genome_extension`  | fna     | File extension of the genome FASTA files.                                               |

A GT-Pro-based abundance estimate (k-mer hits per million reads) is always computed at no extra cost. With `--run_coverm`, `summary/abundance_comparison.tsv` shows per species how well the two agree.

### Genotyping and QC

| Parameter               | Default | Description                                                                                          |
| ----------------------- | ------- | ---------------------------------------------------------------------------------------------------- |
| `--species_bins`        | 64      | Number of parallel species bins (species ID modulo bins).                                            |
| `--pos_offset`          | 1       | Added to GT-Pro global positions in SNV IDs (`<species>_<position>`), as in the paper's catalogue.  |
| `--min_depth`           | 1       | Minimum reads at a site to compute an allele frequency.                                              |
| `--qp_min_depth`        | 2       | Minimum reads at a site to use it for quasi-phaseability.                                            |
| `--qp_site_major_freq`  | 0.9     | A site is "pure" if the major allele has more than this fraction of the reads.                       |
| `--qp_min_frac`         | 0.8     | A sample is quasi-phaseable for a species if more than this fraction of its sites are pure.          |
| `--qp_min_sites`        | 100     | Minimum number of sites to assess quasi-phaseability.                                                |
| `--snv_max_missing`     | 0.9     | Maximum fraction of the cohort's samples without data at an SNV.                                     |
| `--snv_min_maf`         | 0.05    | Minimum minor allele frequency (mean allele frequency over samples with data).                      |
| `--species_min_snvs`    | 1000    | Minimum number of passing SNVs for a species in a cohort.                                            |
| `--min_species_samples` | 10      | Minimum number of samples with the species in a cohort.                                              |

Missingness and MAF are computed per cohort, over all samples of the cohort, as in the paper. An SNV is kept if it passes in at least one cohort in which the species passes. The per-cohort pass flags are in `<species>.snvs.tsv.gz`, so stricter per-cohort selections can be made downstream.

`qp_min_depth`, `qp_min_sites` and `min_species_samples` are not specified in the paper. With `qp_min_depth = 1` every site looks pure, so the default is 2.

### Kinship

| Parameter              | Default | Description                                                         |
| ---------------------- | ------- | ------------------------------------------------------------------- |
| `--skip_kinship`       | false   | Skip the kinship step.                                              |
| `--kinship_qp_only`    | true    | Use only quasi-phaseable samples.                                   |
| `--kinship_min_shared` | 1000    | Minimum SNVs shared by a sample pair (as in the paper).             |
| `--kinship_max_snvs`   | 0       | Subsample SNVs for the distance computation (0 = use all).         |
| `--kinship_n_pcs`      | 10      | Number of principal coordinates written.                            |

## Running on an HPC

Institutional config profiles are loaded automatically, e.g. `-profile snellius` once available. Resources are set per label (`process_single`, `process_low`, `process_medium`, `process_high`) in `conf/base.config` and can be tuned in a custom config with `-c`.

GT-Pro is not on Bioconda/BioContainers. `GTPRO_GENOTYPE` and `GTPRO_SPLITDICT` use the Docker Hub image [`barbarahelena/gt-pro:1.0.1`](https://hub.docker.com/r/barbarahelena/gt-pro) (amd64 and arm64), built from a Dockerfile in a separate container repository. With Singularity/Apptainer, Nextflow converts the image once; set `NXF_SINGULARITY_CACHEDIR` to reuse it across runs.

## Testing

```bash
# wiring test with placeholder data (stub run)
nextflow run . -profile test -stub --outdir results_test
nf-test test tests/default.nf.test
```
