#!/usr/bin/env Rscript
#
# Nextflow template for SUMMARISE.
#
# Outputs
#   sample_summary.tsv          per sample: read pairs, GT-Pro sites, species present / quasi-phaseable
#   species_cohort_summary.tsv  per species x cohort: samples present, quasi-phaseable, SNVs passing QC,
#                               species pass flag, samples retained in the kinship matrix
#   sample_species.tsv.gz       per sample x species coverage, GT-Pro abundance and quasi-phaseability (long)
#   species_abundance_gtpro.tsv sample x species GT-Pro k-mer hits per million reads (abundance proxy)
#   species_abundance.tsv       sample x species relative abundance from CoverM (if run)
#   abundance_comparison.tsv    per species: Spearman correlation of CoverM and GT-Pro abundance (if both)
#   *_mqc.tsv                   MultiQC custom-content tables

suppressPackageStartupMessages(library(data.table))

# Read a (gzipped) table; .gz files are streamed through gzip so that R.utils is not needed
read_table <- function(file, ...) {
    if (endsWith(file, ".gz")) fread(cmd = paste("gzip -dc", shQuote(file)), ...) else fread(file, ...)
}

genome2species_file <- "${genome2species ?: ''}"
samples <- fread("${samples}", colClasses = list(character = c("sample", "cohort")))

read_all <- function(dir, pattern) {
    files <- list.files(dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) return(data.table())
    rbindlist(lapply(files, function(f) {
        d <- read_table(f)
        if (nrow(d) > 0) d[, source_file := basename(f)]
        d
    }), fill = TRUE)
}
strip_suffix <- function(x, suffix) substr(x, 1, nchar(x) - nchar(suffix))

## ---- read counts ---------------------------------------------------------------------------
reads <- read_all("read_stats", "[.]seqkit_stats[.]tsv")
if (nrow(reads) > 0) {
    reads[, sample := strip_suffix(source_file, ".seqkit_stats.tsv")]
    reads <- reads[, .(n_read_files = .N, reads_total = sum(as.numeric(num_seqs)),
                       bases_total = sum(as.numeric(sum_len))), by = sample]
    reads[, read_pairs := reads_total / n_read_files]
}

## ---- per sample x species --------------------------------------------------------------
# GT-Pro abundance proxy: k-mer hits of a species per million reads entering GT-Pro.
# Comparable between samples within a species (not between species: SNV density differs).
sample_species <- read_all("sample_species", "[.]tsv[.]gz")
if (nrow(sample_species) > 0) {
    sample_species[, source_file := NULL]
    if (nrow(reads) > 0 && "gtpro_hits" %in% names(sample_species)) {
        sample_species[reads, reads_total := i.reads_total, on = "sample"]
        sample_species[, gtpro_hits_per_million_reads := 1e6 * gtpro_hits / reads_total]
        sample_species[, reads_total := NULL]
        gt_wide <- dcast(sample_species, sample ~ species, value.var = "gtpro_hits_per_million_reads", fill = 0)
        fwrite(gt_wide, "species_abundance_gtpro.tsv", sep = "\\t")
    }
}
fwrite(sample_species, "sample_species.tsv.gz", sep = "\\t", na = "NA")

## ---- per species x cohort --------------------------------------------------------------
species_cohort <- read_all("species_cohort", "[.]tsv")
if (nrow(species_cohort) > 0) species_cohort[, source_file := NULL]

kin <- read_all("kinship", "[.]kinship_samples[.]tsv")
if (nrow(kin) > 0 && nrow(species_cohort) > 0) {
    kin[, species := as.integer(strip_suffix(source_file, ".kinship_samples.tsv"))]
    kin[samples, cohort := i.cohort, on = "sample"]
    kin_n <- kin[status == "kept", .(n_kinship = .N), by = .(species, cohort)]
    species_cohort <- merge(species_cohort, kin_n, by = c("species", "cohort"), all.x = TRUE)
    species_cohort[is.na(n_kinship), n_kinship := 0L]
}
if (nrow(species_cohort) > 0) setorder(species_cohort, -species_pass, -n_quasi_phaseable, species, cohort)
fwrite(species_cohort, "species_cohort_summary.tsv", sep = "\\t", na = "NA")

## ---- per sample --------------------------------------------------------------------------
summ <- copy(samples)

if (nrow(reads) > 0) summ <- merge(summ, reads, by = "sample", all.x = TRUE)

gtpro <- read_all("gtpro_stats", "[.]gtpro_stats[.]tsv")
if (nrow(gtpro) > 0) {
    gtpro[, source_file := NULL]
    setnames(gtpro, c("n_sites", "n_species", "total_count"),
                    c("gtpro_sites", "gtpro_species", "gtpro_kmer_hits"))
    gtpro[, sample := as.character(sample)]
    summ <- merge(summ, gtpro, by = "sample", all.x = TRUE)
}

if (nrow(sample_species) > 0) {
    sp_n <- sample_species[, .(n_species_present = sum(n_sites > 0),
                               n_species_quasi_phaseable = sum(quasi_phaseable)), by = sample]
    summ <- merge(summ, sp_n, by = "sample", all.x = TRUE)
}
fwrite(summ, "sample_summary.tsv", sep = "\\t", na = "NA")

## ---- species abundance (CoverM) ------------------------------------------------------------
abund_files <- list.files("abundance", pattern = "[.]coverm[.]tsv", full.names = TRUE)
if (length(abund_files) > 0) {
    ab <- rbindlist(lapply(abund_files, function(f) {
        d <- fread(f)
        setnames(d, c("genome", "relative_abundance"))
        d[, sample := strip_suffix(basename(f), ".coverm.tsv")]
        d
    }))
    if (nchar(genome2species_file) > 0) {
        g2s <- read_table(genome2species_file, header = FALSE, col.names = c("genome", "species"),
                     colClasses = "character")
        ab[g2s, species := i.species, on = "genome"]
        ab[is.na(species), species := genome]
    } else {
        ab[, species := genome]
    }
    ab_wide <- dcast(ab, sample ~ species, value.var = "relative_abundance", fun.aggregate = sum, fill = 0)
    fwrite(ab_wide, "species_abundance.tsv", sep = "\\t")

    # agreement between CoverM and the GT-Pro proxy, per species, over samples where both are > 0
    if (nrow(sample_species) > 0 && "gtpro_hits_per_million_reads" %in% names(sample_species)) {
        cm <- ab[, .(relative_abundance = sum(relative_abundance)), by = .(sample, species)]
        gt <- sample_species[, .(sample, species = as.character(species), gtpro_hits_per_million_reads)]
        both <- merge(cm, gt, by = c("sample", "species"))
        cmp <- both[relative_abundance > 0 & gtpro_hits_per_million_reads > 0,
                    .(n_samples = .N,
                      spearman_rho = if (.N >= 10) cor(relative_abundance, gtpro_hits_per_million_reads,
                                                       method = "spearman") else NA_real_),
                    by = species]
        setorder(cmp, -n_samples)
        fwrite(cmp, "abundance_comparison.tsv", sep = "\\t", na = "NA")
    }
}

## ---- MultiQC custom content --------------------------------------------------------------------
mqc_cols <- intersect(c("sample", "cohort", "read_pairs", "gtpro_sites", "gtpro_species",
                        "n_species_present", "n_species_quasi_phaseable"), names(summ))
mqc <- summ[, ..mqc_cols]
writeLines(c(
    "# id: 'gtpro_samples'",
    "# section_name: 'GT-Pro metagenotyping'",
    "# description: 'Reads entering GT-Pro, SNV sites detected, and species present / quasi-phaseable per sample.'",
    "# plot_type: 'table'"
), "gtpro_samples_mqc.tsv")
fwrite(mqc, "gtpro_samples_mqc.tsv", sep = "\\t", append = TRUE, col.names = TRUE, na = "NA")

if (nrow(species_cohort) > 0) {
    top <- species_cohort[species_pass == TRUE]
    if (nrow(top) > 0) {
        top <- top[, .(species_cohort = paste(species, cohort, sep = " | "), n_present,
                       n_quasi_phaseable, n_snvs_pass)]
        writeLines(c(
            "# id: 'species_feasibility'",
            "# section_name: 'Species passing QC'",
            "# description: 'Species x cohort combinations passing SNV QC: samples present, quasi-phaseable samples and SNVs passing QC.'",
            "# plot_type: 'table'"
        ), "species_feasibility_mqc.tsv")
        fwrite(top, "species_feasibility_mqc.tsv", sep = "\\t", append = TRUE, col.names = TRUE)
    }
}

## ---- software versions ----------------------------------------------------------------
writeLines(sprintf("%s\\t%s\\t%s", "${task.process}", c("r-base", "r-data.table"),
                   c(as.character(getRversion()), as.character(packageVersion("data.table")))), "versions.tsv")
