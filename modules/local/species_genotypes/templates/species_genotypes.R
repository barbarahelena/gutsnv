#!/usr/bin/env Rscript
#
# Nextflow template for SPECIES_GENOTYPES (one species bin, all samples).
#
# Build per-species SNV allele-frequency matrices from GT-Pro allele counts,
# assess per-sample strain quasi-phaseability, and apply per-cohort SNV QC.
#
# Inputs
#   counts/*.tsv.gz : GT-Pro allele counts, no header:
#                     species, global_pos, sample, ref_count, alt_count
#   dict            : GT-Pro SNP dictionary rows for this bin, no header:
#                     species, global_pos, ref_allele, alt_allele
#   samples         : sample table with columns sample, cohort (+ subject, timepoint)
#
# Outputs, for every species that passes QC in at least one cohort:
#   <species>.af.tsv.gz       SNV x sample alternative-allele frequency (NA = not observed)
#   <species>.snvs.tsv.gz     SNV annotation + per-cohort missingness / MAF / pass flags
#   <species>.samples.tsv.gz  per-sample coverage, GT-Pro k-mer hits and quasi-phaseability
# and, for all species in the bin (including failing ones):
#   <bin>.sample_species.tsv.gz  per sample x species statistics
#   <bin>.species_cohort.tsv     per species x cohort statistics
#
# Definitions (following Wang et al., Nature 2026, SNV-MWAS):
#   allele frequency  AF = alt_count / (ref_count + alt_count) at sites with depth >= min_depth
#   quasi-phaseable   > qp_min_frac of the sites with depth >= qp_min_depth carry the same
#                     allele in > qp_site_major_freq of reads, over >= qp_min_sites such sites
#   SNV missingness   1 - (samples with data) / (all samples of the cohort)
#   SNV MAF           min(mean AF, 1 - mean AF) over samples with data
#   species pass      >= species_min_snvs SNVs pass in that cohort, and the species is
#                     present in >= min_species_samples samples of that cohort

suppressPackageStartupMessages(library(data.table))

# Read a (gzipped) table; .gz files are streamed through gzip so that R.utils is not needed
read_table <- function(file, ...) {
    if (endsWith(file, ".gz")) fread(cmd = paste("gzip -dc", shQuote(file)), ...) else fread(file, ...)
}

## ---- parameters (interpolated by Nextflow) -----------------------------------------
prefix              <- "bin${bin}"
dict_file           <- "${dict}"
samples_file        <- "${samples}"
pos_offset          <- ${opts.pos_offset}
min_depth           <- ${opts.min_depth}
qp_min_depth        <- ${opts.qp_min_depth}
qp_site_major_freq  <- ${opts.qp_site_major_freq}
qp_min_frac         <- ${opts.qp_min_frac}
qp_min_sites        <- ${opts.qp_min_sites}
min_species_samples <- ${opts.min_species_samples}
snv_max_missing     <- ${opts.snv_max_missing}
snv_min_maf         <- ${opts.snv_min_maf}
species_min_snvs    <- ${opts.species_min_snvs}
setDTthreads(${task.cpus})

## ---- inputs ----------------------------------------------------------------------------
samples <- fread(samples_file, colClasses = list(character = c("sample", "cohort")))
samples[is.na(cohort) | cohort == "", cohort := "all"]
cohort_n <- samples[, .(n_samples_cohort = .N), by = cohort]

count_files <- list.files("counts", pattern = "[.]tsv[.]gz", full.names = TRUE)
if (length(count_files) == 0) stop("No count files found in counts/")

# Stream all per-sample files through a single decompression pipe
file_list <- tempfile()
writeLines(count_files, file_list)
counts <- fread(
    cmd        = sprintf("xargs gzip -dc < %s", shQuote(file_list)),
    header     = FALSE,
    col.names  = c("species", "global_pos", "sample", "ref_count", "alt_count"),
    colClasses = c("integer", "integer", "character", "integer", "integer")
)
counts <- counts[sample %in% samples[["sample"]]]
counts[, depth := ref_count + alt_count]
# total GT-Pro k-mer hits per species and sample (abundance proxy), before any depth filter
hits <- counts[, .(gtpro_hits = sum(depth)), keyby = .(species, sample)]
counts <- counts[depth >= min(min_depth, qp_min_depth)]
setkey(counts, species)

dict <- read_table(
    dict_file, header = FALSE,
    col.names  = c("species", "global_pos", "ref_allele", "alt_allele"),
    colClasses = c("integer", "integer", "character", "character")
)
setkey(dict, species, global_pos)

## ---- per species ---------------------------------------------------------------------
sample_species_out <- list()
species_cohort_out <- list()

for (sp in sort(unique(counts[["species"]]))) {

    x <- counts[.(sp)]

    ## per-sample coverage and quasi-phaseability
    cov <- x[depth >= min_depth, .(n_sites = .N, mean_depth = mean(depth)), by = sample]
    qp  <- x[depth >= qp_min_depth,
             .(n_qp_sites = .N,
               frac_major = mean(pmax(ref_count, alt_count) / depth > qp_site_major_freq)),
             by = sample]
    ss <- merge(cov, qp, by = "sample", all = TRUE)
    ss <- merge(ss, hits[.(sp), .(sample, gtpro_hits)], by = "sample", all = TRUE)
    ss[is.na(n_sites), n_sites := 0L]
    ss[is.na(gtpro_hits), gtpro_hits := 0L]
    ss[is.na(n_qp_sites), n_qp_sites := 0L]
    ss[, quasi_phaseable := !is.na(frac_major) & n_qp_sites >= qp_min_sites & frac_major > qp_min_frac]
    ss <- merge(samples[, .(sample, cohort)], ss, by = "sample")
    ss[, species := sp]
    setcolorder(ss, c("sample", "cohort", "species"))
    sample_species_out[[length(sample_species_out) + 1]] <- ss

    ## per-cohort SNV QC on all samples with data (paper definition)
    g <- x[depth >= min_depth, .(global_pos, sample, af = alt_count / depth)]
    g[samples, cohort := i.cohort, on = "sample"]

    snv_stats <- g[, .(n_obs = .N, af_mean = mean(af)), by = .(global_pos, cohort)]
    snv_stats[cohort_n, n_samples_cohort := i.n_samples_cohort, on = "cohort"]
    snv_stats[, missing := 1 - n_obs / n_samples_cohort]
    snv_stats[, maf := pmin(af_mean, 1 - af_mean)]
    snv_stats[, snv_pass := missing <= snv_max_missing & maf >= snv_min_maf]

    sc <- ss[, .(n_present = sum(n_sites > 0), n_quasi_phaseable = sum(quasi_phaseable)), by = cohort]
    sc <- merge(cohort_n, sc, by = "cohort", all.x = TRUE)
    sc[is.na(n_present), `:=`(n_present = 0L, n_quasi_phaseable = 0L)]
    sc <- merge(sc,
                snv_stats[, .(n_snvs_observed = .N, n_snvs_pass = sum(snv_pass)), by = cohort],
                by = "cohort", all.x = TRUE)
    sc[is.na(n_snvs_observed), `:=`(n_snvs_observed = 0L, n_snvs_pass = 0L)]
    sc[, species_pass := n_snvs_pass >= species_min_snvs & n_present >= min_species_samples]
    sc[, species := sp]
    setcolorder(sc, c("species", "cohort"))
    species_cohort_out[[length(species_cohort_out) + 1]] <- sc

    pass_cohorts <- sc[species_pass == TRUE, cohort]
    if (length(pass_cohorts) == 0) next

    # an SNV is kept if it passes in at least one cohort in which the species passes
    snv_stats[, snv_pass := snv_pass & cohort %in% pass_cohorts]
    keep_pos <- snv_stats[snv_pass == TRUE, unique(global_pos)]
    g <- g[global_pos %in% keep_pos]

    ## SNV annotation table
    snvs <- g[, .(n_obs = .N, af_mean = mean(af)), by = global_pos]
    snvs[, maf := pmin(af_mean, 1 - af_mean)]
    wide_stats <- dcast(
        snv_stats[global_pos %in% keep_pos],
        global_pos ~ cohort,
        value.var = c("n_obs", "missing", "maf", "snv_pass")
    )
    snvs <- merge(snvs, wide_stats, by = "global_pos", all.x = TRUE)
    snvs <- merge(dict[.(sp), .(global_pos, ref_allele, alt_allele)], snvs, by = "global_pos", all.y = TRUE)
    snvs[, `:=`(species = sp, pos = global_pos + pos_offset)]
    snvs[, snv_id := paste(species, pos, sep = "_")]
    setcolorder(snvs, c("snv_id", "species", "pos", "global_pos", "ref_allele", "alt_allele"))
    setorder(snvs, pos)

    ## SNV x sample allele-frequency matrix
    g[, af := round(af, 4)]
    af <- dcast(g, global_pos ~ sample, value.var = "af")
    af <- af[match(snvs[["global_pos"]], global_pos)]
    af[, global_pos := NULL]
    af <- cbind(snv_id = snvs[["snv_id"]], af)

    fwrite(af,   sprintf("%s.af.tsv.gz", sp),      sep = "\\t", na = "NA")
    fwrite(snvs, sprintf("%s.snvs.tsv.gz", sp),    sep = "\\t", na = "NA")
    fwrite(ss,   sprintf("%s.samples.tsv.gz", sp), sep = "\\t", na = "NA")
}

## ---- bin-level summaries (all species) -----------------------------------------------
fwrite(rbindlist(sample_species_out, fill = TRUE),
       sprintf("%s.sample_species.tsv.gz", prefix), sep = "\\t", na = "NA")
fwrite(rbindlist(species_cohort_out, fill = TRUE),
       sprintf("%s.species_cohort.tsv", prefix), sep = "\\t", na = "NA")

## ---- software versions ----------------------------------------------------------------
writeLines(sprintf("%s\\t%s\\t%s", "${task.process}", c("r-base", "r-data.table"),
                   c(as.character(getRversion()), as.character(packageVersion("data.table")))), "versions.tsv")
