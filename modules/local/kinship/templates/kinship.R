#!/usr/bin/env Rscript
#
# Nextflow template for KINSHIP (one species).
#
# Strain-level genetic distance, kinship (relatedness) matrix and PCoA,
# following Wang et al., Nature 2026 (SNV-MWAS):
#
#   D_ij = mean |AF_i - AF_j| over SNVs observed in both samples
#          (equivalent to the paper's Manhattan distance on dosages, scaled by 2 x shared SNVs)
#   pairs sharing < min_shared SNVs are set to missing, and the sample with the most
#   missing pairs is removed repeatedly until no missing values remain
#   K_ij = 1 - D_ij
#   PCoA on D with Cailliez correction (ape::pcoa)
#
# Inputs : af       SNV x sample allele-frequency matrix (first column snv_id)
#          samples  per-sample quasi-phaseability for this species
# Outputs: <species>.kinship.tsv.gz, <species>.distance.tsv.gz, <species>.nshared.tsv.gz,
#          <species>.pcs.tsv, <species>.pcoa_eigenvalues.tsv, <species>.kinship_samples.tsv

suppressPackageStartupMessages({
    library(data.table)
    library(ape)
})

# Read a (gzipped) table; .gz files are streamed through gzip so that R.utils is not needed
read_table <- function(file, ...) {
    if (endsWith(file, ".gz")) fread(cmd = paste("gzip -dc", shQuote(file)), ...) else fread(file, ...)
}

## ---- parameters (interpolated by Nextflow) -----------------------------------------
prefix     <- "${species}"
af_file    <- "${af}"
ss_file    <- "${samples}"
qp_only    <- ${opts.qp_only ? 'TRUE' : 'FALSE'}
min_shared <- ${opts.min_shared}
max_snvs   <- ${opts.max_snvs}
n_pcs      <- ${opts.n_pcs}
block_size <- 200L
setDTthreads(${task.cpus})
set.seed(1)

## ---- software versions ----------------------------------------------------------------
writeLines(sprintf("%s\\t%s\\t%s", "${task.process}", c("r-base", "r-data.table", "r-ape"),
                   c(as.character(getRversion()), as.character(packageVersion("data.table")), as.character(packageVersion("ape")))), "versions.tsv")

## ---- inputs ----------------------------------------------------------------------------
af <- read_table(af_file, na.strings = "NA")
n_snvs_total <- nrow(af)
X <- as.matrix(af[, -1])
rm(af)

ss <- read_table(ss_file)
log_tab <- data.table(sample = colnames(X), status = "kept")
log_tab[ss, quasi_phaseable := i.quasi_phaseable, on = "sample"]

if (qp_only) {
    keep <- colnames(X) %in% ss[quasi_phaseable == TRUE, sample]
    log_tab[!keep, status := "not_quasi_phaseable"]
    X <- X[, keep, drop = FALSE]
}

if (max_snvs > 0 && nrow(X) > max_snvs) {
    X <- X[sort(sample.int(nrow(X), max_snvs)), , drop = FALSE]
}
log_tab[, n_snvs_used := nrow(X)]

finish_without_kinship <- function(msg) {
    message(msg)
    fwrite(log_tab, sprintf("%s.kinship_samples.tsv", prefix), sep = "\\t")
    quit(save = "no", status = 0)
}
if (ncol(X) < 3) finish_without_kinship("Fewer than 3 samples available; no kinship computed.")

## ---- shared SNVs and iterative sample removal ------------------------------------------
observed <- !is.na(X)
nshared  <- crossprod(observed + 0)
dimnames(nshared) <- list(colnames(X), colnames(X))

ok <- nshared >= min_shared
repeat {
    n_bad <- colSums(!ok)
    if (all(n_bad == 0)) break
    worst <- names(which.max(n_bad))
    log_tab[sample == worst, status := "too_few_shared_snvs"]
    ok <- ok[rownames(ok) != worst, colnames(ok) != worst, drop = FALSE]
}
kept <- rownames(ok)
if (length(kept) < 3) finish_without_kinship("Fewer than 3 samples left after the shared-SNV filter.")

X       <- X[, kept, drop = FALSE]
nshared <- nshared[kept, kept]

## ---- distance: mean absolute AF difference over shared SNVs ----------------------------
n <- ncol(X)
D <- matrix(0, n, n, dimnames = list(kept, kept))
for (i in seq_len(n - 1)) {
    j_all <- (i + 1):n
    for (start in seq(1, length(j_all), by = block_size)) {
        j <- j_all[start:min(start + block_size - 1, length(j_all))]
        s <- colSums(abs(X[, j, drop = FALSE] - X[, i]), na.rm = TRUE)
        D[i, j] <- s / nshared[i, j]
        D[j, i] <- D[i, j]
    }
}
K <- 1 - D

## ---- PCoA --------------------------------------------------------------------------------
pc      <- pcoa(D, correction = "cailliez")
use_cor <- !is.null(pc[["vectors.cor"]])
vec     <- if (use_cor) pc[["vectors.cor"]] else pc[["vectors"]]
rel     <- if (use_cor) pc[["values"]][["Rel_corr_eig"]] else pc[["values"]][["Relative_eig"]]
k       <- min(n_pcs, ncol(vec))
pcs     <- data.table(sample = kept, vec[, seq_len(k), drop = FALSE])
setnames(pcs, c("sample", paste0("PC", seq_len(k))))
eig     <- data.table(pc = paste0("PC", seq_along(rel)), relative_eigenvalue = rel)

## ---- outputs -----------------------------------------------------------------------------
write_mat <- function(m, file, digits = 6) {
    fwrite(data.table(sample = rownames(m), round(m, digits)), file, sep = "\\t")
}
write_mat(K, sprintf("%s.kinship.tsv.gz", prefix))
write_mat(D, sprintf("%s.distance.tsv.gz", prefix))
write_mat(nshared, sprintf("%s.nshared.tsv.gz", prefix), digits = 0)
fwrite(pcs, sprintf("%s.pcs.tsv", prefix), sep = "\\t")
fwrite(eig[seq_len(min(nrow(eig), 50))], sprintf("%s.pcoa_eigenvalues.tsv", prefix), sep = "\\t")
fwrite(log_tab, sprintf("%s.kinship_samples.tsv", prefix), sep = "\\t")
