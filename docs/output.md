# gutsnv: Output

All outputs are published by the workflow `output` block in `main.nf` into `--outdir`.

## Genotypes (`genotypes/`)

One directory per species that passes QC in at least one cohort. `genotypes/index.csv` lists the files per species, for downstream scripts.

| File                       | Content                                                                                                                                                                  |
| -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `<species>.af.tsv.gz`      | SNV × sample alternative-allele frequency (0–1, `NA` = no reads). First column `snv_id`.                                                                                 |
| `<species>.snvs.tsv.gz`    | Per SNV: `snv_id` (`<species>_<global_pos + pos_offset>`), GT-Pro `global_pos`, `ref_allele`, `alt_allele`, overall `n_obs`/`af_mean`/`maf`, and per cohort `n_obs_*`, `missing_*`, `maf_*`, `snv_pass_*`. |
| `<species>.samples.tsv.gz` | Per sample: cohort, `n_sites`, `mean_depth`, `n_qp_sites`, `frac_major` (fraction of pure sites), `gtpro_hits` (k-mer hits) and `quasi_phaseable`.                  |

For association analyses, use only quasi-phaseable samples of a species. The allele frequency of such samples is close to 0 or 1.

### Species abundance as a covariate

Two abundance measures are available, both with exactly the same species definition as the genotypes (the GT-Pro reference genomes):

- **CoverM** (`species_abundance.tsv`): relative abundance from read mapping to one reference genome per species, as in Wang et al. 2026. Optional, with `--run_coverm`.
- **GT-Pro proxy** (`species_abundance_gtpro.tsv`): k-mer hits of a species per million reads. Always computed, at no extra cost (the default abundance measure). It scales with the species' coverage, so it is comparable between samples *within* a species (which is what a per-species covariate needs), but not between species, because the number of SNVs in the catalogue differs per species.

`abundance_comparison.tsv` shows how well the two agree per species.

In the wide tables the column names are numeric species IDs; read them with an explicit header (e.g. `data.table::fread(file, header = TRUE)`).

## Kinship (`kinship/`)

One directory per species, with index `kinship/index.csv`.

| File                                | Content                                                                     |
| ----------------------------------- | --------------------------------------------------------------------------- |
| `<species>.kinship.tsv.gz`          | Sample × sample kinship, K = 1 − D.                                         |
| `<species>.distance.tsv.gz`         | Sample × sample distance D (mean absolute allele-frequency difference).     |
| `<species>.nshared.tsv.gz`          | SNVs shared by each sample pair.                                            |
| `<species>.pcs.tsv`                 | Principal coordinates (Cailliez-corrected PCoA of D).                       |
| `<species>.pcoa_eigenvalues.tsv`    | Relative eigenvalues.                                                       |
| `<species>.kinship_samples.tsv`     | Per sample: `kept`, `not_quasi_phaseable` or `too_few_shared_snvs`.         |

The kinship matrix is written as computed (as in the paper). It is not forced to be positive semi-definite.

## Summary (`summary/`)

| File                         | Content                                                                                                                       |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `species_cohort_summary.tsv` | Per species × cohort: samples present, quasi-phaseable, SNVs observed/passing, `species_pass`, samples kept in the kinship.  |
| `sample_summary.tsv`         | Per sample: read pairs, GT-Pro sites and species, species present and quasi-phaseable.                                        |
| `sample_species.tsv.gz`      | Per sample × species coverage, GT-Pro abundance and quasi-phaseability, for all species (also those failing QC).               |
| `species_abundance_gtpro.tsv` | Sample × species GT-Pro k-mer hits per million reads: an abundance proxy with exactly the genotyped species definition.   |
| `species_abundance.tsv`      | Sample × species relative abundance (CoverM), if run.                                                                         |
| `abundance_comparison.tsv`   | Per species: Spearman correlation between CoverM and the GT-Pro proxy (samples where both > 0; ≥ 10 samples), if CoverM ran. |

## QC and other outputs

| Directory          | Content                                                                    |
| ------------------ | -------------------------------------------------------------------------- |
| `qc/fastqc/`       | FastQC reports of the raw reads.                                           |
| `qc/fastp/`        | fastp reports.                                                             |
| `qc/host_removal/` | Bowtie2 logs (host alignment rate).                                        |
| `qc/read_stats/`   | seqkit statistics of the reads entering GT-Pro.                            |
| `gtpro/stats/`     | Per sample: GT-Pro sites, species and k-mer hits.                          |
| `gtpro/raw/`       | Raw GT-Pro output (with `--save_gtpro_raw`).                               |
| `abundance/`       | CoverM output per sample (with `--run_coverm`).                            |
| `clean_reads/`     | Reads entering genotyping (with `--save_clean_reads`).                     |
| `multiqc/`         | MultiQC report and `software_versions.yml`.                                |
| `pipeline_info/`   | Nextflow execution report, timeline, trace and DAG.                        |
