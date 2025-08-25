#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Prepare files for calculating C+T based PRS
# Author: Shelley
# Date: 2025-08-03
# Usage: Rscript 03_prepare_C+T_files.R <CV_number>
# ------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(RISCA)
  library(boot)
})

# ---- Input parameter ----
args <- commandArgs(trailingOnly = TRUE)
CV_number <- as.character(args[1])
cat(sprintf("Processing CV %s at %s.\n", CV_number, Sys.time()))

# ---- Directories (edit BASE_DIR as needed) ----
BASE_DIR <- "PATH/TO/PROJECT"

gwas_dir   <- file.path(BASE_DIR, "CV-GWAS/Step2")
clump_dir  <- file.path(BASE_DIR, "cvPGS/PRS_methods/C+T/clumps")
out_dir    <- file.path(BASE_DIR, "cvPGS/PRS_methods/C+T/C+T_PRS", CV_number)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Read GWAS summary statistics for chr1–22 ----
chr_vec <- 1:22
gwas_files <- sprintf("%s/regenie_step2_CV%s_chr%d_TS_log_Zadj.regenie", gwas_dir, CV_number, chr_vec)

cat("Reading GWAS summary statistics for all chromosomes...\n")
gwas_list <- lapply(seq_along(gwas_files), function(i) {
  cat(sprintf("  Reading chromosome %2d: %s\n", chr_vec[i], gwas_files[i]))
  fread(gwas_files[i], data.table = FALSE)
})
cat("Finished reading GWAS summary statistics.\n")

dat <- bind_rows(gwas_list)
cat(sprintf("Combined GWAS data: %d variants from all chromosomes.\n", nrow(dat)))

colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","INFO","N",
                   "TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
dat$P <- 10^(-1*dat$LOG10P)
dat <- dat[,c("CHROM","ID","REF","POS","ALT","BETA","P")]
colnames(dat) <- c("CHR","ID","REF","BP","ALT","BETA","P")

# ---- 2. Read clumped files for chr1–22 ----
clumped_files <- sprintf("%s/LD_clump_CV%s_chr%d.clumped", clump_dir, CV_number, chr_vec)

cat("Reading LD clumped files for all chromosomes...\n")
clump_list <- lapply(seq_along(clumped_files), function(i) {
  cat(sprintf("  Reading clumped file for chr %2d: %s\n", chr_vec[i], clumped_files[i]))
  fread(clumped_files[i], data.table = FALSE)
})
cat("Finished reading all clumped files.\n")

LD <- bind_rows(clump_list)
cat(sprintf("Combined LD clumped data: %d variants from all chromosomes.\n", nrow(LD)))

# ---- 3. Extract index SNPs and join with summary data ----
clump.snp <- LD[,3,drop=FALSE]
dat$SNP <- paste(dat$CHR, dat$BP, dat$REF, dat$ALT, sep="_")

prs.all <- left_join(clump.snp, dat, by = c("SNP" = "ID"))
cat(sprintf("Number of index SNPs in clump.snp: %d\n", nrow(clump.snp)))
cat(sprintf("Number of index SNPs in prs.all: %d\n", nrow(prs.all)))

if (nrow(prs.all) != nrow(clump.snp)) {
  stop("Mismatch in number of index SNPs between clump.snp and prs.all.")
}

# ---- 4. Write PRS coefficient & p-value files ----
prs.file <- prs.all[, c("SNP", "ALT", "BETA")]
p.value.file <- prs.all[, c("SNP", "P")]

write.table(prs.file,  file.path(out_dir, "C+T_prs_coeff.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(p.value.file, file.path(out_dir, "C+T_p_value.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

beta_final <- inner_join(prs.file, p.value.file, by="SNP")
write.table(beta_final, file.path(out_dir, "C+T_prs_coeff_and_pvalue.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

cat(sprintf("Number of variants in beta_final: %d\n", nrow(beta_final)))

# ---- 5. Save p-value thresholds (q_range) ----
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)

q_range <- data.frame(
  name = paste0("p_value_", seq_along(pthres)),
  from = 0,
  to   = pthres
)

write.table(q_range, file.path(out_dir, "C+T_q_range.txt"),
            sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

cat("Preview of q_range:\n"); print(q_range)
cat(sprintf("Saved q_range to: %s\n", file.path(out_dir, "C+T_q_range.txt")))

# ---- Done ----
cat(sprintf("C+T PRS files for CV %s have been generated at %s.\n", 
            CV_number, Sys.time()))
