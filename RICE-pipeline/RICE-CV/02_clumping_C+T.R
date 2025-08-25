#!/usr/bin/env Rscript
# This script is used to clump GWAS summary statistics for C+T methods
# Author: Shelley
# Date: 2025-07-27
# Usage: Rscript 02_clumping_C+T.R <CV_number> <chr>

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript 02_clumping_C+T.R <CV_number> <chr>")
}
CV_number <- args[1]   # e.g., "01"
chr <- args[2]         # e.g., "7"

# ---- File paths (edit BASE_DIR as needed) ----
BASE_DIR <- "PATH/TO/PROJECT"

regenie_output <- sprintf("%s/regenie_step2_CV%s_chr%s_TS_log_Zadj.regenie", 
                          file.path(BASE_DIR, "CV-GWAS/Step2"), CV_number, chr)

assoc_out <- sprintf("%s/association.results_CV%s_chr%s_assoc.txt", 
                     file.path(BASE_DIR, "cvPGS/PRS_methods/C+T"), CV_number, chr)

pheno_model_file <- file.path(BASE_DIR, "Sample.Splitting/pheno_model.tsv")
keep_file        <- file.path(BASE_DIR, "Sample.Splitting/pheno_model.keep")

plink   <- "plink"   # ensure PLINK is in your PATH
bfile   <- sprintf("%s/chr%s.imputation.reset.varid", 
                   file.path(BASE_DIR, "Imputation/reset_VarID"), chr)

clump_out <- sprintf("%s/LD_clump_CV%s_chr%s", 
                     file.path(BASE_DIR, "cvPGS/PRS_methods/C+T/clumps"), CV_number, chr)

# Clumping thresholds
pthr   <- 1
r2thr  <- 0.1
kbpthr <- 500

# ---- 1. Read REGENIE output and prepare summary stats ----
dat <- read.table(regenie_output, header = TRUE, stringsAsFactors = FALSE)
colnames(dat) <- c("CHROM","POS","ID","REF","ALT","A1_FREQ","INFO","N",
                   "TEST","BETA","SE","CHISQ","LOG10P","EXTRA")
dat$P <- 10^(-1 * dat$LOG10P)

# Unique SNP ID: chr_pos_ref_alt
dat$VARID <- paste0(dat$CHROM, "_", dat$POS, "_", dat$REF, "_", dat$ALT)
dat <- dat[, c("CHROM", "VARID", "REF", "POS", "ALT", "BETA", "P")]
colnames(dat) <- c("CHR", "SNP", "REF", "BP", "ALT", "BETA", "P")

write.table(dat, file = assoc_out, col.names = TRUE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

# ---- 2. Prepare PLINK --keep file (FID IID) ----
if (!file.exists(keep_file)) {
  pheno_model <- read.table(pheno_model_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  write.table(data.frame(FID=pheno_model[,1], IID=pheno_model[,1]),
              file = keep_file, col.names=FALSE, row.names=FALSE, 
              quote=FALSE, sep="\t")
}

# ---- 3. Run PLINK C+T clumping ----
cmd <- sprintf('%s --bfile %s --keep %s --clump %s --clump-p1 %s --clump-r2 %s --clump-kb %s --threads 8 --out %s',
               plink, bfile, keep_file, assoc_out, pthr, r2thr, kbpthr, clump_out)

cat("[INFO] Running PLINK command:\n", cmd, "\n")
system(cmd)