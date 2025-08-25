#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Pipeline to get LDpred2 and LASSOsum2 estimates
# NOTE: You must pre-compute LD matrices using LDstore2.
#       The script assumes that LD matrices are saved as:
#       LD_with_blocks_chr{1..22}.rds
#
# Author: Shelley
# Date: 2025-07-31
# Usage: Rscript 05_LDpred2_LASSOsum2.R <CV_number>
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(pROC)
  library(bigsnpr)
  library(bigsparser)
  library(readr)
  library(boot)
  library(bigstatsr)
  library(RISCA)
  library(bigreadr)
})

# ---- Input ----
args <- commandArgs(trailingOnly = TRUE)
CV_number <- as.character(args[1])

# ---- Paths ----
BASE_DIR   <- "PATH/TO/PROJECT"
map_ldref  <- file.path(BASE_DIR, "LD_block_from_LDpred2/map_hm3_plus.rds")
gwas_dir   <- file.path(BASE_DIR, "CV-GWAS/Step2")
ld_dir     <- file.path(BASE_DIR, "LD_block_from_LDpred2")   # contains LD_with_blocks_chr*.rds
tmpdir     <- file.path(BASE_DIR, "LDpred2_LASSOsum2/tmp-data", CV_number)
outdir     <- file.path(BASE_DIR, "LDpred2_LASSOsum2/LDpred2_LASSOsum2_estimates", CV_number)
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Read reference map ----
map_ldref <- readRDS(map_ldref)
cat("Reading HapMap3+ variants...\n")
cat("Number of variants in map file:", nrow(map_ldref), "\n")

# ---- 2. Read GWAS summary statistics (chr1–22) ----
gwas_list <- lapply(1:22, function(chr) {
  path <- sprintf("%s/regenie_step2_CV%s_chr%d_TS_log_Zadj.regenie", gwas_dir, CV_number, chr)
  cat(sprintf("[CV %s] Reading chromosome %2d: %s\n", CV_number, chr, path))
  if (!file.exists(path)) stop("File not found: ", path)
  fread2(path)
})
gwas_all <- bind_rows(gwas_list)
cat("Total SNPs in GWAS read:", nrow(gwas_all), "\n")

# Reformat
colnames(gwas_all) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","INFO","N","TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
gwas_all$P <- 10^(-1*gwas_all$LOG10P)
sumstats <- gwas_all[,c('CHROM','ID','POS','REF','ALT','BETA','SE','P','N',"A1_FREQ","INFO")]
names(sumstats) <- c("chr","rsid","pos","a0","a1","beta","beta_se","p","n_eff","A1_FREQ","INFO")

# ---- 3. Match SNPs with HapMap3+ reference ----
info_snp <- snp_match(sumstats, map_ldref)
rownames(info_snp) <- info_snp$rsid
cat("Number of variants after matching:", nrow(info_snp), "\n")

# ---- 4. SNP QC ----
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 1 / sqrt(n_eff * beta_se^2 + beta^2))
is_bad <- sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.05 | sd_ldref < 0.05
df_beta <- info_snp[!is_bad, ]
cat("SNPs remaining after QC:", nrow(df_beta), "\n")

# ---- 5. Build genome-wide correlation matrix (LD matrix) ----
# NOTE: LDstore2 must be used beforehand to generate LD_with_blocks_chr*.rds
tmp <- tempfile(tmpdir = tmpdir)
for (chr in 1:22) {
  cat("Adding LD for chr", chr, "...\n")
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  corr_chr <- readRDS(file.path(ld_dir, sprintf("LD_with_blocks_chr%d.rds", chr)))[ind.chr3, ind.chr3]
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}
cat("Successfully built genome-wide correlation matrix!\n")

# ---- 6. LD score regression ----
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                               chi2 = (beta / beta_se)^2,
                               sample_size = n_eff,
                               ncores = 1))
print(ldsc)
h2_est <- ldsc[["h2"]]

# ---- 7. LDpred2 ----
h2_seq <- seq(0.1, 1.5, by = 0.1)
p_seq  <- signif(seq_log(1e-4, 1, length.out = 17), 2)
params <- expand.grid(p = p_seq, h2 = signif(abs(h2_est) * h2_seq, 3), sparse = c(FALSE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, burn_in = 50, num_iter = 200, ncores = 1)
rownames(beta_grid) <- df_beta$rsid
beta_grid <- cbind(beta_grid, df_beta[,c('chr','pos','a0','a1','rsid')])
beta_grid[is.na(beta_grid)] <- 0
beta_grid <- as.data.frame(beta_grid)
cat("LDpred2 completed.\n")

# ---- 8. LASSOsum2 ----
delta_path <- function(max=100, min=0.5, n=10){
  sqrt_max <- max^(1/3); sqrt_min <- min^(1/3)
  path <- numeric(n)
  for (i in 1:n) {
    path[n+1-i] = (sqrt_max - (sqrt_max-sqrt_min)/(n-1)  * (i-1) )^3;
  }
  return(path)
}
beta_LASSOsum2 <- snp_lassosum2(corr, df_beta, delta = delta_path(), ncores = 1, maxiter=1000)
params2 <- attr(beta_LASSOsum2, "grid_param")
rownames(beta_LASSOsum2) <- df_beta$rsid
beta_LASSOsum2 <- cbind(beta_LASSOsum2, df_beta[,c('chr','pos','a0','a1','rsid')])
beta_LASSOsum2[is.na(beta_LASSOsum2)] <- 0
beta_LASSOsum2 <- as.data.frame(beta_LASSOsum2)
cat("LASSOsum2 completed.\n")

# ---- 9. Save outputs ----
saveRDS(params, file = file.path(outdir, sprintf("params_LDpred2_CV%s.rds", CV_number)))
saveRDS(beta_grid, file = file.path(outdir, sprintf("beta_grid_LDpred2_CV%s.rds", CV_number)))
saveRDS(params2, file = file.path(outdir, sprintf("params_LASSOsum2_CV%s.rds", CV_number)))
saveRDS(beta_LASSOsum2, file = file.path(outdir, sprintf("beta_LASSOsum2_CV%s.rds", CV_number)))

prs.file  <- data.frame(SNP = beta_grid$rsid,  CHR=beta_grid$chr, POS=beta_grid$pos, REF=beta_grid$a0, ALT=beta_grid$a1, BETA = beta_grid[,1:(ncol(beta_grid)-5)])
prs.file2 <- data.frame(SNP = beta_LASSOsum2$rsid, CHR=beta_LASSOsum2$chr, POS=beta_LASSOsum2$pos, REF=beta_LASSOsum2$a0, ALT=beta_LASSOsum2$a1, BETA = beta_LASSOsum2[,1:(ncol(beta_LASSOsum2)-5)])

write.table(prs.file,  file = file.path(outdir, sprintf("LDpred2_CV%s_prs.txt", CV_number)),  col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(prs.file2, file = file.path(outdir, sprintf("LASSOsum2_CV%s_prs.txt", CV_number)), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

cat("All LDpred2 and LASSOsum2 estimates saved for CV", CV_number, "\n")
