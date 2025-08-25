#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# Phenotype processing + 10-fold CV splits + REGENIE train files
# Author: Shelley | Date: 2025-07-26 (v3, simplified)
# ------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

# ===================== EDIT THESE PATHS ===========================
WES_FAM    <- "PATH/TO/chr1.fam"                 # WES .fam (QCed)
IMP_SAMPLE <- "PATH/TO/ukb22828_c1_b0_v3.sample" # Imputation .sample (WTCHG)
PHENO_FILE <- "PATH/TO/participant_all.tsv"      # UKB phenotype TSV
OUTDIR     <- "PATH/TO/output_dir"               # Output directory
# ================================================================

# Basic settings
SEED      <- 2025
KFOLDS    <- 10
PCS_N     <- 10                       # number of PCs to keep (p22009_a1..aN)
PHENO_COL <- "p22192_i0"              # TL: Z-adjusted log T/S

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n")

# ---------------- Module 1: shared genotype sample IDs ----------------
msg("Reading WES .fam: %s", WES_FAM)
wes_fam <- fread(WES_FAM, header = FALSE, data.table = FALSE)
wes_ids <- wes_fam[, 2]

msg("Reading imputation .sample: %s", IMP_SAMPLE)
imp_sample <- fread(IMP_SAMPLE, header = FALSE, skip = 2, data.table = FALSE)
imp_ids <- imp_sample[, 1]

intersect_ids <- intersect(wes_ids, imp_ids)
msg("Intersect(WES, Imputation) = %d samples", length(intersect_ids))

# ---------------- Module 2: phenotype prep ---------------------------
msg("Reading phenotype table: %s", PHENO_FILE)
pheno <- fread(PHENO_FILE)
stopifnot("eid" %in% names(pheno))

pheno_shared <- pheno[eid %in% intersect_ids]
msg("Participants with genotype & phenotype rows: %d", nrow(pheno_shared))

# Define phenotype & covariates
pcs_cols   <- paste0("p22009_a", 1:PCS_N)
covar_cols <- c("p21022", "p31", pcs_cols, "p21000_i0")  # age, sex, PCs, self-report ancestry
need_cols  <- c("eid", PHENO_COL, covar_cols)

missing <- setdiff(need_cols, names(pheno_shared))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

pheno_model <- pheno_shared[, ..need_cols]
setnames(pheno_model, "p21022", "age")
pheno_model[, age2 := age^2]

# Friendly names
setnames(
  pheno_model,
  old = c("eid", PHENO_COL, "age", "age2", "p31", "p21000_i0", pcs_cols),
  new = c("eid", "TS_log_Zadj", "age", "age2", "sex", "self_report_ancestry",
          paste0("PC", 1:PCS_N))
)

na_n <- pheno_model[is.na(TS_log_Zadj), .N]
msg("NA in TS_log_Zadj: %d", na_n)
pheno_model <- pheno_model[!is.na(TS_log_Zadj)]
msg("Samples after removing NA TL: %d", nrow(pheno_model))

# ---------------- Module 3: 10-fold CV splits -----------------------
set.seed(SEED)
n <- nrow(pheno_model)

folds <- sample(rep(1:KFOLDS, length.out = n))
pheno_model[, fold := folds]

all_val_ids <- integer(0)
for (k in 1:KFOLDS) {
  val_ids  <- pheno_model[fold == k, eid]
  tune_ids <- pheno_model[fold == ifelse(k < KFOLDS, k + 1, 1), eid]
  train_ids <- pheno_model[!(eid %in% c(val_ids, tune_ids)), eid]

  stopifnot(length(intersect(val_ids, tune_ids)) == 0,
            length(intersect(val_ids, train_ids)) == 0,
            length(intersect(tune_ids, train_ids)) == 0)

  fwrite(data.table(eid = train_ids), file.path(OUTDIR, sprintf("cv%02d_train.txt", k)), col.names = FALSE)
  fwrite(data.table(eid = tune_ids),  file.path(OUTDIR, sprintf("cv%02d_tune.txt",  k)), col.names = FALSE)
  fwrite(data.table(eid = val_ids),   file.path(OUTDIR, sprintf("cv%02d_val.txt",   k)), col.names = FALSE)

  all_val_ids <- c(all_val_ids, val_ids)
  msg("CV %02d: train=%d  tune=%d  val=%d", k, length(train_ids), length(tune_ids), length(val_ids))
}

stopifnot(length(all_val_ids) == n,
          length(unique(all_val_ids)) == n,
          setequal(all_val_ids, pheno_model$eid))

# Save the processed phenotype table
fwrite(pheno_model, file.path(OUTDIR, "pheno_model.tsv"), sep = "\t", quote = FALSE)
msg("Saved: %s", file.path(OUTDIR, "pheno_model.tsv"))

# ---------------- Module 4: REGENIE training files -------------------
covar_keep <- c("age", "age2", "sex", paste0("PC", 1:PCS_N))
stopifnot(all(covar_keep %in% names(pheno_model)))

for (k in 1:KFOLDS) {
  train_ids <- fread(file.path(OUTDIR, sprintf("cv%02d_train.txt", k)), header = FALSE)[[1]]
  train_dt  <- pheno_model[eid %in% train_ids]

  # Phenotype
  pheno_out <- train_dt[, .(FID = eid, IID = eid, TS_log_Zadj)]
  fwrite(pheno_out, file.path(OUTDIR, sprintf("cv%02d_train.pheno.txt", k)),
         sep = "\t", col.names = TRUE, na = "NA", quote = FALSE)

  # Covariates
  covar_out <- train_dt[, c(.(FID = eid, IID = eid), .SD), .SDcols = covar_keep]
  fwrite(covar_out, file.path(OUTDIR, sprintf("cv%02d_train.covariate.txt", k)),
         sep = "\t", col.names = TRUE, na = "NA", quote = FALSE)
}

msg("Done. Outputs written to: %s", OUTDIR)
# ------------------------------------------------------------------