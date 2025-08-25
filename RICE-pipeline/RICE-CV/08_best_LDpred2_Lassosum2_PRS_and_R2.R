#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Rscript to get the best LDPred2 and LASSOSum2 PRS and R2 in the tuning dataset
# Author: Shelley
# Date: 2025-08-02
# Usage: Rscript 08_best_LDpred2_Lassosum2_PRS_and_R2.R <CV_number>
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
})

# ---------------- Args ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide the CV number as the first argument.")
CV_number <- as.character(args[1])

# ---------------- Base dirs (EDIT THIS) ----------------
BASE_DIR <- "PATH/TO/PROJECT"

LDpred2_dir   <- file.path(BASE_DIR, "LDpred2_LASSOsum2/LDpred2_LASSOsum2_estimates/PRS/LDpred2",   CV_number)
LASSOsum2_dir <- file.path(BASE_DIR, "LDpred2_LASSOsum2/LDpred2_LASSOsum2_estimates/PRS/LASSOsum2", CV_number)

# Per-chromosome PRS .sscore files
chr_vec <- 1:22
prs_train_files  <- sprintf("%s/CV%s_chr%d_prs_train.sscore",      LDpred2_dir,   CV_number, chr_vec)
prs_tune_files   <- sprintf("%s/CV%s_chr%d_prs_tune.sscore",       LDpred2_dir,   CV_number, chr_vec)
prs_vali_files   <- sprintf("%s/CV%s_chr%d_prs_validation.sscore", LDpred2_dir,   CV_number, chr_vec)

lasso_train_files <- sprintf("%s/CV%s_chr%d_prs_train.sscore",      LASSOsum2_dir, CV_number, chr_vec)
lasso_tune_files  <- sprintf("%s/CV%s_chr%d_prs_tune.sscore",       LASSOsum2_dir, CV_number, chr_vec)
lasso_vali_files  <- sprintf("%s/CV%s_chr%d_prs_validation.sscore", LASSOsum2_dir, CV_number, chr_vec)

# ---------------- Read LDpred2 PRS .sscore ----------------
cat("Reading LDpred2 PRS files...\n")
prs_train_list <- lapply(seq_along(prs_train_files), function(i) {
  cat(sprintf("  [CV %s] LDpred2 train chr %2d\n", CV_number, chr_vec[i])); fread(prs_train_files[i], data.table = FALSE)
})
prs_tune_list  <- lapply(seq_along(prs_tune_files), function(i) {
  cat(sprintf("  [CV %s] LDpred2 tune  chr %2d\n", CV_number, chr_vec[i])); fread(prs_tune_files[i],  data.table = FALSE)
})
prs_vali_list  <- lapply(seq_along(prs_vali_files), function(i) {
  cat(sprintf("  [CV %s] LDpred2 vali  chr %2d\n", CV_number, chr_vec[i])); fread(prs_vali_files[i],  data.table = FALSE)
})

# ---------------- Read LASSOsum2 PRS .sscore ----------------
cat("Reading LASSOsum2 PRS files...\n")
lasso_train_list <- lapply(seq_along(lasso_train_files), function(i) {
  cat(sprintf("  [CV %s] LASSOsum2 train chr %2d\n", CV_number, chr_vec[i])); fread(lasso_train_files[i], data.table = FALSE)
})
lasso_tune_list  <- lapply(seq_along(lasso_tune_files), function(i) {
  cat(sprintf("  [CV %s] LASSOsum2 tune  chr %2d\n", CV_number, chr_vec[i])); fread(lasso_tune_files[i],  data.table = FALSE)
})
lasso_vali_list  <- lapply(seq_along(lasso_vali_files), function(i) {
  cat(sprintf("  [CV %s] LASSOsum2 vali  chr %2d\n", CV_number, chr_vec[i])); fread(lasso_vali_files[i],  data.table = FALSE)
})

# ---------------- Helpers ----------------
add_chr_suffix <- function(df, chr) {
  score_cols <- grep("^SCORE[0-9]+_SUM$", colnames(df), value = TRUE)
  colnames(df)[match(score_cols, colnames(df))] <- paste0(score_cols, "_chr", chr)
  df
}

# LDpred2 expects 255 scores (SCORE1..SCORE255)
get_genomewide_prs_matrix_LDpred2 <- function(prs_list, chr_vec, label = "LDpred2") {
  labeled_list <- lapply(seq_along(prs_list), function(i) add_chr_suffix(prs_list[[i]], chr_vec[i]))
  prs_all <- Reduce(function(x, y) merge(x, y, by = c("#FID", "IID")), labeled_list)
  score_cols_chr1 <- grep("^SCORE[0-9]+_SUM_chr1$", colnames(prs_all), value = TRUE)
  base_scores <- sub("_chr1$", "", score_cols_chr1)
  if (length(base_scores) != 255)
    stop(sprintf("Expected 255 SCORE_SUM columns for %s, found %d", label, length(base_scores)))
  prs_mat <- prs_all[, c("#FID", "IID")]
  for (score in base_scores) {
    cols <- paste0(score, "_chr", chr_vec)
    prs_mat[[score]] <- rowSums(prs_all[, cols, drop = FALSE])
  }
  rownames(prs_mat) <- prs_mat$IID
  prs_mat
}

# LASSOsum2 expects 300 scores (SCORE1..SCORE300)
get_genomewide_prs_matrix_LASSOsum2 <- function(prs_list, chr_vec, label = "LASSOsum2") {
  labeled_list <- lapply(seq_along(prs_list), function(i) add_chr_suffix(prs_list[[i]], chr_vec[i]))
  prs_all <- Reduce(function(x, y) merge(x, y, by = c("#FID", "IID")), labeled_list)
  score_cols_chr1 <- grep("^SCORE[0-9]+_SUM_chr1$", colnames(prs_all), value = TRUE)
  base_scores <- sub("_chr1$", "", score_cols_chr1)
  if (length(base_scores) != 300)
    stop(sprintf("Expected 300 SCORE_SUM columns for %s, found %d", label, length(base_scores)))
  prs_mat <- prs_all[, c("#FID", "IID")]
  for (score in base_scores) {
    cols <- paste0(score, "_chr", chr_vec)
    prs_mat[[score]] <- rowSums(prs_all[, cols, drop = FALSE])
  }
  rownames(prs_mat) <- prs_mat$IID
  prs_mat
}

# ---------------- Build genome-wide PRS matrices ----------------
LDpred2_train_prs_mat <- get_genomewide_prs_matrix_LDpred2(prs_train_list, chr_vec, "LDpred2 train")
LDpred2_tune_prs_mat  <- get_genomewide_prs_matrix_LDpred2(prs_tune_list,  chr_vec, "LDpred2 tune")
LDpred2_vali_prs_mat  <- get_genomewide_prs_matrix_LDpred2(prs_vali_list,  chr_vec, "LDpred2 validation")

LASSOsum2_train_prs_mat <- get_genomewide_prs_matrix_LASSOsum2(lasso_train_list, chr_vec, "LASSOsum2 train")
LASSOsum2_tune_prs_mat  <- get_genomewide_prs_matrix_LASSOsum2(lasso_tune_list,  chr_vec, "LASSOsum2 tune")
LASSOsum2_vali_prs_mat  <- get_genomewide_prs_matrix_LASSOsum2(lasso_vali_list,  chr_vec, "LASSOsum2 validation")

cat("LDpred2 train/tune/vali dims:", dim(LDpred2_train_prs_mat), dim(LDpred2_tune_prs_mat), dim(LDpred2_vali_prs_mat), "\n")
cat("LASSOsum2 train/tune/vali dims:", dim(LASSOsum2_train_prs_mat), dim(LASSOsum2_tune_prs_mat), dim(LASSOsum2_vali_prs_mat), "\n")

# ---------------- Save PRS matrices ----------------
dir.create(LDpred2_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(LASSOsum2_dir, recursive = TRUE, showWarnings = FALSE)

save(LDpred2_train_prs_mat, file = file.path(LDpred2_dir,   "LDpred2_train_prs_mat.RData"))
save(LDpred2_tune_prs_mat,  file = file.path(LDpred2_dir,   "LDpred2_tune_prs_mat.RData"))
save(LDpred2_vali_prs_mat,  file = file.path(LDpred2_dir,   "LDpred2_vali_prs_mat.RData"))

save(LASSOsum2_train_prs_mat, file = file.path(LASSOsum2_dir, "LASSOsum2_train_prs_mat.RData"))
save(LASSOsum2_tune_prs_mat,  file = file.path(LASSOsum2_dir, "LASSOsum2_tune_prs_mat.RData"))
save(LASSOsum2_vali_prs_mat,  file = file.path(LASSOsum2_dir, "LASSOsum2_vali_prs_mat.RData"))

cat("Saved PRS matrices as .RData files.\n")

# ---------------- Phenotypes & splits ----------------
split_dir  <- file.path(BASE_DIR, "Sample.Splitting")
pheno_file <- file.path(split_dir, "pheno_model.tsv")

pheno_all <- fread(pheno_file, data.table = FALSE)
train_ids <- fread(file.path(split_dir, sprintf("cv%s_train.txt", CV_number)), header = FALSE)[[1]]
tune_ids  <- fread(file.path(split_dir, sprintf("cv%s_tune.txt",  CV_number)), header = FALSE)[[1]]
vali_ids  <- fread(file.path(split_dir, sprintf("cv%s_val.txt",   CV_number)), header = FALSE)[[1]]

pheno_train <- pheno_all %>% filter(eid %in% train_ids) %>% mutate(IID = eid)
pheno_tune  <- pheno_all %>% filter(eid %in% tune_ids)  %>% mutate(IID = eid)
pheno_vali  <- pheno_all %>% filter(eid %in% vali_ids)  %>% mutate(IID = eid)

LDpred2_train_data <- left_join(pheno_train, LDpred2_train_prs_mat, by = "IID")
LDpred2_tune_data  <- left_join(pheno_tune,  LDpred2_tune_prs_mat,  by = "IID")
LDpred2_vali_data  <- left_join(pheno_vali,  LDpred2_vali_prs_mat,  by = "IID")

LASSOsum2_train_data <- left_join(pheno_train, LASSOsum2_train_prs_mat, by = "IID")
LASSOsum2_tune_data  <- left_join(pheno_tune,  LASSOsum2_tune_prs_mat,  by = "IID")
LASSOsum2_vali_data  <- left_join(pheno_vali,  LASSOsum2_vali_prs_mat,  by = "IID")

cat(sprintf("n(train/tune/vali) after merge (LDpred2): %d / %d / %d\n",
            nrow(LDpred2_train_data), nrow(LDpred2_tune_data), nrow(LDpred2_vali_data)))
cat(sprintf("n(train/tune/vali) after merge (LASSOsum2): %d / %d / %d\n",
            nrow(LASSOsum2_train_data), nrow(LASSOsum2_tune_data), nrow(LASSOsum2_vali_data)))

# ---------------- R² for LDpred2 (255 scores) ----------------
h2_seq <- seq(0.1, 1.5, by = 0.1)
p_seq  <- signif(seq_log(1e-4, 1, length.out = 17), 2)
sets_LDpred2 <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE))

r2_tun_vec_LDpred2 <- rep(0, nrow(sets_LDpred2))
model.null <- lm(TS_log_Zadj ~ age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                 data = LDpred2_tune_data)
for (k in seq_len(nrow(sets_LDpred2))) {
  prs <- LDpred2_tune_data[!is.na(LDpred2_tune_data$TS_log_Zadj), paste0("SCORE", k, "_SUM")]
  model.prs <- lm(model.null$residuals ~ prs, data = LDpred2_tune_data)
  r2_tun_vec_LDpred2[k] <- summary(model.prs)$r.squared
}
idx_LDpred2 <- which.max(r2_tun_vec_LDpred2)

best_prs_train_LDpred2 <- LDpred2_train_data[, c("IID", paste0("SCORE", idx_LDpred2, "_SUM"))]
best_prs_tune_LDpred2  <- LDpred2_tune_data[,  c("IID", paste0("SCORE", idx_LDpred2, "_SUM"))]
best_prs_validation_LDpred2 <- LDpred2_vali_data[, c("IID", paste0("SCORE", idx_LDpred2, "_SUM"))]

# ---------------- R² for LASSOsum2 (300 scores) ----------------
r2_tun_vec_LASSOsum2 <- rep(0, 300)
model.null <- lm(TS_log_Zadj ~ age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                 data = LASSOsum2_tune_data)
for (k in 1:300) {
  prs <- LASSOsum2_tune_data[!is.na(LASSOsum2_tune_data$TS_log_Zadj), paste0("SCORE", k, "_SUM")]
  model.prs <- lm(model.null$residuals ~ prs, data = LASSOsum2_tune_data)
  r2_tun_vec_LASSOsum2[k] <- summary(model.prs)$r.squared
}
idx_LASSOsum2 <- which.max(r2_tun_vec_LASSOsum2)

best_prs_train_LASSOsum2 <- LASSOsum2_train_data[, c("IID", paste0("SCORE", idx_LASSOsum2, "_SUM"))]
best_prs_tune_LASSOsum2  <- LASSOsum2_tune_data[,  c("IID", paste0("SCORE", idx_LASSOsum2, "_SUM"))]
best_prs_validation_LASSOsum2 <- LASSOsum2_vali_data[, c("IID", paste0("SCORE", idx_LASSOsum2, "_SUM"))]

cat(sprintf("Best R2 (LDpred2, tuning): %.4f (SCORE%d_SUM)\n", max(r2_tun_vec_LDpred2), idx_LDpred2))
cat(sprintf("Best R2 (LASSOsum2, tuning): %.4f (SCORE%d_SUM)\n", max(r2_tun_vec_LASSOsum2), idx_LASSOsum2))

# ---------------- Save best PRS, R2, and indices ----------------
dir.create(LDpred2_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(LASSOsum2_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(best_prs_train_LDpred2,      file.path(LDpred2_dir,   "best_prs_train_LDpred2.txt"),      sep = "\t")
fwrite(best_prs_tune_LDpred2,       file.path(LDpred2_dir,   "best_prs_tune_LDpred2.txt"),       sep = "\t")
fwrite(best_prs_validation_LDpred2, file.path(LDpred2_dir,   "best_prs_validation_LDpred2.txt"), sep = "\t")

fwrite(best_prs_train_LASSOsum2,      file.path(LASSOsum2_dir, "best_prs_train_LASSOsum2.txt"),      sep = "\t")
fwrite(best_prs_tune_LASSOsum2,       file.path(LASSOsum2_dir, "best_prs_tune_LASSOsum2.txt"),       sep = "\t")
fwrite(best_prs_validation_LASSOsum2, file.path(LASSOsum2_dir, "best_prs_validation_LASSOsum2.txt"), sep = "\t")

saveRDS(r2_tun_vec_LDpred2,   file.path(LDpred2_dir,   "r2_tune_vec_LDpred2.rds"))
saveRDS(r2_tun_vec_LASSOsum2, file.path(LASSOsum2_dir, "r2_tune_vec_LASSOsum2.rds"))

fwrite(data.frame(idx_LDpred2 = idx_LDpred2),       file.path(LDpred2_dir,   "best_prs_index_LDpred2.txt"),
       col.names = FALSE, row.names = FALSE, quote = FALSE)
fwrite(data.frame(idx_LASSOsum2 = idx_LASSOsum2),   file.path(LASSOsum2_dir, "best_prs_index_LASSOsum2.txt"),
       col.names = FALSE, row.names = FALSE, quote = FALSE)

cat(sprintf("Best PRS and R2 saved. Tuning complete for CV %s at %s.\n", CV_number, Sys.time()))
