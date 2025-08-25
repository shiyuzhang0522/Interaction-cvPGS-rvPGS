#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Rscript to get the optimal C+T PRS in the tuning dataset
# Author: Shelley
# Date: 2025-08-04
# Usage: Rscript 07_CT_optimal_PRS.R <CV_number>
# ------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(RISCA)
  library(boot)
})

# ---------------- Args ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide the CV number as the first argument.")
CV_number <- as.character(args[1])
cat(sprintf("Processing for CV %s at %s.\n", CV_number, Sys.time()))

# ---------------- Parameters ----------------
pthres <- c(5E-08,5E-07,5E-06,5E-05,5E-04,5E-03,5E-02,5E-01,1.0)
n_pthres <- length(pthres)
chr_vec <- 1:22

# ---------------- Base dirs (edit as needed) ----------------
BASE_DIR <- "PATH/TO/PROJECT"
prs_dir  <- file.path(BASE_DIR, "cvPGS/PRS_methods/C+T/C+T_PRS", CV_number)
split_dir <- file.path(BASE_DIR, "Sample.Splitting")
pheno_file <- file.path(split_dir, "pheno_model.tsv")

# ---------------- Read .sscore files ----------------
read_sscore_list <- function(set, CV_number, n_pthres, chr_vec) {
  prs_per_p <- vector("list", n_pthres)
  file_count <- 0
  for (pthres_idx in seq_len(n_pthres)) {
    prs_per_p[[pthres_idx]] <- vector("list", length(chr_vec))
    for (chr in chr_vec) {
      sscore_file <- sprintf("%s/CV%s.chr%d.C+T_prs_%s.p_value_%d.sscore",
                             prs_dir, CV_number, chr, set, pthres_idx)
      prs_per_p[[pthres_idx]][[chr]] <- fread(sscore_file)
      file_count <- file_count + 1
      cat(sprintf("Read: %s\n", sscore_file))
    }
  }
  expected <- n_pthres * length(chr_vec)
  cat(sprintf("Sanity check: Read %d %s files, expected %d.\n",
              file_count, set, expected))
  if (file_count != expected) stop("File count mismatch for ", set, " set.")
  prs_per_p
}

train_prs_per_chr_list <- read_sscore_list("train", CV_number, n_pthres, chr_vec)
tune_prs_per_chr_list  <- read_sscore_list("tune", CV_number, n_pthres, chr_vec)
vali_prs_per_chr_list  <- read_sscore_list("validation", CV_number, n_pthres, chr_vec)

# ---------------- Helper to merge genome-wide PRS ----------------
add_chr_suffix <- function(df, chr) {
  score_col <- grep("^SCORE[0-9]+_SUM$", colnames(df), value = TRUE)
  df <- as.data.frame(df[, c("#FID", "IID", score_col), with = FALSE])
  colnames(df)[3] <- paste0(score_col, "_chr", chr)
  df
}
sum_genomewide_prs <- function(per_chr_list, chr_vec) {
  labeled_list <- lapply(seq_along(per_chr_list), function(i)
    add_chr_suffix(per_chr_list[[i]], chr_vec[i]))
  merged <- Reduce(function(x, y) merge(x, y, by = c("#FID", "IID")), labeled_list)
  score_cols <- grep("^SCORE[0-9]+_SUM_chr[0-9]+$", colnames(merged), value = TRUE)
  merged[score_cols] <- lapply(merged[score_cols], as.numeric)
  prs_mat <- merged[, c("#FID", "IID")]
  prs_mat$PRS <- rowSums(merged[, score_cols, drop = FALSE], na.rm = TRUE)
  prs_mat
}

# ---------------- Collapse per p-threshold ----------------
collapse_sets <- function(per_chr_list, ids, setname) {
  prs_cols <- vector("list", n_pthres)
  for (i in seq_len(n_pthres)) {
    prs <- sum_genomewide_prs(per_chr_list[[i]], chr_vec)
    prs_cols[[i]] <- prs[, "PRS", drop = FALSE]
    colnames(prs_cols[[i]]) <- paste0("p_value_", i)
  }
  out <- cbind(per_chr_list[[1]][[1]][, c("#FID", "IID")], do.call(cbind, prs_cols))
  fwrite(out, file.path(prs_dir, sprintf("C+T_prs_all_%s.txt", setname)), sep = "\t")
  out
}

prs_mat_train <- collapse_sets(train_prs_per_chr_list, chr_vec, "train")
prs_mat_tune  <- collapse_sets(tune_prs_per_chr_list,  chr_vec, "tune")
prs_mat_vali  <- collapse_sets(vali_prs_per_chr_list,  chr_vec, "validation")

# ---------------- Phenotype & covariates ----------------
pheno_all <- fread(pheno_file, data.table = FALSE)
train_ids <- fread(file.path(split_dir, sprintf("cv%s_train.txt", CV_number)), header = FALSE)[[1]]
tune_ids  <- fread(file.path(split_dir, sprintf("cv%s_tune.txt",  CV_number)), header = FALSE)[[1]]
vali_ids  <- fread(file.path(split_dir, sprintf("cv%s_val.txt",   CV_number)), header = FALSE)[[1]]

pheno_train <- pheno_all %>% filter(eid %in% train_ids) %>% mutate(IID = eid)
pheno_tune  <- pheno_all %>% filter(eid %in% tune_ids)  %>% mutate(IID = eid)
pheno_vali  <- pheno_all %>% filter(eid %in% vali_ids)  %>% mutate(IID = eid)

pheno_train <- left_join(pheno_train, prs_mat_train, by = "IID")
pheno_tune  <- left_join(pheno_tune,  prs_mat_tune,  by = "IID")
pheno_vali  <- left_join(pheno_vali,  prs_mat_vali,  by = "IID")

# ---------------- R² per p-threshold ----------------
r2_tun_vec <- rep(0, length(pthres))
model.null <- lm(TS_log_Zadj ~ age + age2 + sex + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                 data = pheno_tune)
for (k in 1:length(pthres)) {
  prs <- pheno_tune[!is.na(pheno_tune$TS_log_Zadj), paste0("p_value_", k)]
  model.prs <- lm(model.null$residuals ~ prs, data = pheno_tune)
  r2_tun_vec[k] <- summary(model.prs)$r.squared
}

idx_CT <- which.max(r2_tun_vec)
cat(sprintf("Best R² for CV %s: %.6f at p-value threshold: %g (index: %d)\n",
            CV_number, r2_tun_vec[idx_CT], pthres[idx_CT], idx_CT))

# ---------------- Save best PRS ----------------
save_best <- function(pheno_set, setname, idx_CT) {
  df <- pheno_set[, c("IID", "#FID", paste0("p_value_", idx_CT))]
  colnames(df) <- c("IID", "FID", "PRS")
  fwrite(df, file = file.path(prs_dir, sprintf("C+T_prs_best_%s.txt", setname)), sep = "\t")
}

save_best(pheno_train, "train", idx_CT)
save_best(pheno_tune,  "tune",  idx_CT)
save_best(pheno_vali,  "validation", idx_CT)

cat(sprintf("Best C+T PRS for CV %s has been saved at %s.\n", CV_number, Sys.time()))
