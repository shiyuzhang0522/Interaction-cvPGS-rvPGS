#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Rscript to get the ensemble rvPGS (keep all significant masks across 10 folds)
# Author: Shelley
# Date: 2025-08-07
# Usage: Rscript 02_ensemble_rvPGS.R <CV_number>
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(gdsfmt)
  library(SeqArray)
  library(SeqVarTools)
  library(GENESIS)
  library(STAAR)
  library(STAARpipeline)
  library(STAARpipelineSummary)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(data.table)
  library(dplyr)
  library(RISCA)
  library(boot)
  library(stringr)
  library(caret)
  library(ranger)
  library(glmnet)
})

# ---------------- Args ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Please provide the CV number as the first argument.")
CV_number <- as.character(args[1])
cat(sprintf("Calculating ensemble rvPGS for CV %s at %s.\n", CV_number, Sys.time()))

# ---------------- Base directory (EDIT THIS) ----------------
BASE_DIR <- "PATH/TO/PROJECT"

# ---------------- Paths (filenames preserved) ----------------
# Function to build G* burden per gene/mask
gstar_fn      <- file.path(BASE_DIR, "RICE-pipeline/RICE-RV/01_Gene_Centric_Coding_G_Star.R")

# Significant training p-values (per fold)
train_pvals_file <- file.path(
  BASE_DIR,
  "RICE-RV/rvPGS/coding_results",
  sprintf("CV%s_UKBB_WES_Coding_Train_coding_sig.QCed.csv", CV_number)
)

# agds directory pointer (RData object that evaluates to a named vector indexed by chr)
agds_dir_rdata <- file.path(
  BASE_DIR,
  "RICE-RV/STAARpipeline/pre-steps/gds2agds/agds_dir.Rdata"
)

# Null model per fold (Train)
null_model_path <- file.path(
  BASE_DIR,
  "RICE-RV/STAARpipeline/Step1.Nullmodels",
  CV_number,
  sprintf("CV%s_Train_Null_Model.RData", CV_number)
)

# Annotation catalog (data.frame with columns name, dir)
anno_catalog_rdata <- file.path(
  BASE_DIR,
  "RICE-RV/STAARpipeline/pre-steps/gds2agds/Annotation_name_catalog.Rdata"
)

# Phenotypes / splits
split_dir  <- file.path(BASE_DIR, "Sample.Splitting")
pheno_file <- file.path(split_dir, "pheno_model.tsv")
train_id_file <- file.path(split_dir, sprintf("cv%s_train.txt", CV_number))
tune_id_file  <- file.path(split_dir, sprintf("cv%s_tune.txt",  CV_number))
vali_id_file  <- file.path(split_dir, sprintf("cv%s_val.txt",   CV_number))

# Output
out_dir <- file.path(BASE_DIR, "RICE-RV/rvPGS/Ensemble_rvPGS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- Source burden function ----------------
source(gstar_fn)  # defines Gene_Centric_Coding_G_Star()

# ---------------- Load significant gene/mask list ----------------
Train_PVals_All <- read.csv(train_pvals_file)
Train_PVals_All <- Train_PVals_All[Train_PVals_All$STAARB <= 1e-3, ]
cat(sprintf("Loaded %d gene-masks with STAARB <= 1e-3 for CV%s\n",
            nrow(Train_PVals_All), CV_number))

if (nrow(Train_PVals_All) == 0) {
  stop("No significant gene-masks at STAARB <= 1e-3. Cannot build rvPGS.")
}

# ---------------- Load AGDS dir & null model ----------------
agds_dir <- get(load(agds_dir_rdata))
obj_nullmodel <- get(load(null_model_path))

# ---------------- Phenotypes & splits ----------------
pheno_all <- fread(pheno_file, data.table = FALSE)

train_ids <- fread(train_id_file, header = FALSE)[[1]]
tune_ids  <- fread(tune_id_file,  header = FALSE)[[1]]
vali_ids  <- fread(vali_id_file,  header = FALSE)[[1]]

pheno_train <- pheno_all %>% filter(eid %in% train_ids) %>% mutate(IID = eid)
pheno_tune  <- pheno_all %>% filter(eid %in% tune_ids)  %>% mutate(IID = eid)
pheno_valid <- pheno_all %>% filter(eid %in% vali_ids)  %>% mutate(IID = eid)

# Use all split IDs for extraction
obj_nullmodel$id_include <- c(pheno_train$IID, pheno_tune$IID, pheno_valid$IID)
ids_gstar <- obj_nullmodel$id_include

# ---------------- Parameters for extraction ----------------
QC_label <- "annotation/info/QC_label"
geno_missing_imputation <- "mean"
variant_type <- "SNV"
Annotation_dir <- "annotation/info/FunctionalAnnotation"
Annotation_name_catalog <- get(load(anno_catalog_rdata))

cat(sprintf("Extracting G* genotype burden matrix at %s…\n", Sys.time()))

# ---------------- Build G* (per significant gene/mask) ----------------
G_star_list <- vector("list", nrow(Train_PVals_All))

for (i in seq_len(nrow(Train_PVals_All))) {
  chr       <- Train_PVals_All$Chr[i]
  gene_name <- Train_PVals_All$Gene[i]
  category  <- Train_PVals_All$Category[i]

  gds_path <- agds_dir[as.character(chr)]
  if (is.na(gds_path) || !file.exists(gds_path)) {
    stop("AGDS file not found for chr ", chr, ": ", gds_path)
  }

  genofile <- seqOpen(gds_path)
  G_star_list[[i]] <- Gene_Centric_Coding_G_Star(
    chr              = chr,
    gene_name        = gene_name,
    category         = category,
    genofile         = genofile,
    obj_nullmodel    = obj_nullmodel,
    rare_maf_cutoff  = 0.01,
    rv_num_cutoff    = 2,
    QC_label         = QC_label,
    variant_type     = variant_type,
    geno_missing_imputation = geno_missing_imputation,
    Annotation_dir   = Annotation_dir,
    Annotation_name_catalog = Annotation_name_catalog,
    # pass gene coordinates table explicitly if your function expects it:
    # genes_info = YOUR_GENE_TABLE,
    silent           = TRUE
  )
  seqClose(genofile)
  if (i %% 50 == 0) cat("  processed ", i, " / ", nrow(Train_PVals_All), "\n")
}

# Combine to matrix
G_star_gene_centric_coding <- do.call(cbind, G_star_list)
cat(sprintf("Dimensions of G*: %d rows x %d columns\n",
            nrow(G_star_gene_centric_coding), ncol(G_star_gene_centric_coding)))

# Round burdens (imputation may yield non-integers)
G_star_gene_centric_coding <- round(G_star_gene_centric_coding)

# ---------------- Filter columns:
# (1) >10 samples carry ≥1 rare allele
# (2) Total rare allele count >10
keep_mask <- (apply(G_star_gene_centric_coding, 2, function(x) sum(x != 0)) > 10) &
             (colSums(G_star_gene_centric_coding) > 10)
cat("Columns kept after burden filters: ", sum(keep_mask), " / ", ncol(G_star_gene_centric_coding), "\n")
if (!any(keep_mask)) stop("No burden columns pass filters; cannot proceed.")

G_star_gene_centric_coding <- G_star_gene_centric_coding[, keep_mask, drop = FALSE]
Train_PVals_All <- Train_PVals_All[keep_mask, , drop = FALSE]

# ---------------- Split G* by dataset ----------------
G_train <- G_star_gene_centric_coding[ids_gstar %in% pheno_train$IID, , drop = FALSE]
G_tune  <- G_star_gene_centric_coding[ids_gstar %in% pheno_tune$IID,  , drop = FALSE]
G_val   <- G_star_gene_centric_coding[ids_gstar %in% pheno_valid$IID, , drop = FALSE]

# Build design matrices aligned with phenotypes
build_X <- function(pheno_df, G_mat) {
  dat <- data.frame(IID = pheno_df$IID, G_mat, check.names = FALSE)
  colnames(dat) <- c("IID", paste0("X", seq_len(ncol(G_mat))))
  pheno_df2 <- dplyr::inner_join(pheno_df, dat, by = "IID")
  X <- as.matrix(pheno_df2[, paste0("X", seq_len(ncol(G_mat))), drop = FALSE])
  list(pheno = pheno_df2, X = X)
}
bt <- build_X(pheno_train, G_train); pheno_train <- bt$pheno; X_train <- bt$X
bt <- build_X(pheno_tune,  G_tune);  pheno_tune  <- bt$pheno; X_tune  <- bt$X
bt <- build_X(pheno_valid, G_val);   pheno_validation <- bt$pheno; X_valid <- bt$X

# ---------------- Residualize TL on covariates (null models) ----------------
covars <- ~ age + age2 + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
fml <- as.formula(paste("TS_log_Zadj", deparse(covars)))

y_train      <- lm(fml, data = pheno_train)$residuals
y_tune       <- lm(fml, data = pheno_tune)$residuals
y_validation <- lm(fml, data = pheno_validation)$residuals

pheno_train$y_train <- NA_real_;      pheno_train$y_train[!is.na(pheno_train$TS_log_Zadj)]      <- y_train
pheno_tune$y_tune   <- NA_real_;      pheno_tune$y_tune[!is.na(pheno_tune$TS_log_Zadj)]          <- y_tune
pheno_validation$y_validation <- NA_real_; pheno_validation$y_validation[!is.na(pheno_validation$TS_log_Zadj)] <- y_validation

# ---------------- Train base models (train set) ----------------
trait <- "TS_log_Zadj"
idx_tr <- !is.na(pheno_train[[trait]])

lasso_train <- glmnet(X_train[idx_tr, , drop = FALSE], pheno_train$y_train[idx_tr], family = "gaussian", alpha = 1)
ridge_train <- glmnet(X_train[idx_tr, , drop = FALSE], pheno_train$y_train[idx_tr], family = "gaussian", alpha = 0)
lm_train    <- lm.fit(cbind(1, X_train[idx_tr, , drop = FALSE]), pheno_train$y_train[idx_tr])
lm_train$coefficients[is.na(lm_train$coefficients)] <- 0

# Save per-mask betas for reference
beta_matrix <- as.data.frame(cbind(
  as.matrix(lasso_train$beta),
  as.matrix(ridge_train$beta),
  lm_train$coefficients[-1]
))
colnames(beta_matrix) <- c(
  paste0("lasso_prs", seq_len(ncol(as.matrix(lasso_train$beta)))),
  paste0("ridge_prs", seq_len(ncol(as.matrix(ridge_train$beta)))),
  "lm_prs1"
)
beta_matrix <- cbind(Train_PVals_All[, 1:4, drop = FALSE], beta_matrix)

# ---------------- PRS predictions for tune/validation ----------------
lasso_prs_tune <- as.matrix(predict(lasso_train, X_tune))
ridge_prs_tune <- as.matrix(predict(ridge_train, X_tune))
lm_prs_tune    <- as.numeric(cbind(1, X_tune) %*% matrix(lm_train$coefficients, ncol = 1))

lasso_prs_val  <- as.matrix(predict(lasso_train, X_valid))
ridge_prs_val  <- as.matrix(predict(ridge_train, X_valid))
lm_prs_val     <- as.numeric(cbind(1, X_valid) %*% matrix(lm_train$coefficients, ncol = 1))

all_prs_tune  <- as.data.frame(cbind(lasso_prs_tune, ridge_prs_tune, lm_prs_tune))
colnames(all_prs_tune) <- c(paste0("lasso_prs", ncol(lasso_prs_tune) %>% seq_len()),
                            paste0("ridge_prs", ncol(ridge_prs_tune) %>% seq_len()),
                            "lm_prs1")
all_prs_val   <- as.data.frame(cbind(lasso_prs_val, ridge_prs_val, lm_prs_val))
colnames(all_prs_val)  <- c(paste0("lasso_prs", ncol(lasso_prs_val) %>% seq_len()),
                            paste0("ridge_prs", ncol(ridge_prs_val) %>% seq_len()),
                            "lm_prs1")

# ---------------- Drop unusable / highly correlated columns ----------------
# (1) drop columns that are entirely NA after correlation attempt
mtx <- suppressWarnings(cor(all_prs_tune))
drop <- names(all_prs_tune)[apply(mtx, 2, function(x) sum(is.na(x))) >= (nrow(mtx) - 1)]
if (length(drop)) {
  all_prs_tune <- dplyr::select(all_prs_tune, -all_of(drop))
  all_prs_val  <- dplyr::select(all_prs_val,  -all_of(drop))
}

# (2) drop high-correlation > 0.98
if (ncol(all_prs_tune) > 1) {
  mtx2 <- suppressWarnings(cor(all_prs_tune))
  drop2_idx <- findCorrelation(mtx2, cutoff = 0.98)
  if (length(drop2_idx)) {
    drop2 <- names(all_prs_tune)[drop2_idx]
    all_prs_tune <- dplyr::select(all_prs_tune, -all_of(drop2))
    all_prs_val  <- dplyr::select(all_prs_val,  -all_of(drop2))
  }
}

# (3) drop perfect linear combos
lc <- caret::findLinearCombos(all_prs_tune)
if (!is.null(lc$remove) && length(lc$remove)) {
  drop3 <- names(all_prs_tune)[lc$remove]
  all_prs_tune <- dplyr::select(all_prs_tune, -all_of(drop3))
  all_prs_val  <- dplyr::select(all_prs_val,  -all_of(drop3))
}

# ---------------- Simple ensemble selection on tuning ----------------
# Compare linear fits of each candidate to y_tune; keep best (highest R^2),
# then regress y on that best + original features to get final coefficients,
# matching your original approach.
y <- pheno_tune$y_tune
Xcand <- as.matrix(all_prs_tune[!is.na(y), , drop = FALSE])
y     <- y[!is.na(y)]

if (ncol(Xcand) == 0) stop("No candidate PRS features remain after filtering.")

R2_vec <- apply(Xcand, 2, function(col) {
  summary(lm(y ~ col))$r.squared
})

best_idx <- which.max(R2_vec)
best_col <- Xcand[, best_idx, drop = FALSE]
ens_fit  <- lm(y ~ ., data = data.frame(best_col))
ens_pred <- ens_fit$fitted.values

# Regress y on ensemble prediction + original features to get final weights
final_coefs <- coef(lm(y ~ ., data = data.frame(y = ens_pred, Xcand)))
final_coefs[is.na(final_coefs)] <- 0

# ---------------- Convert base betas to final weights ----------------
# Map final coefficients back to the beta_matrix columns (names must match)
nm <- c("(Intercept)", colnames(all_prs_tune))
names(final_coefs) <- ifelse(is.na(names(final_coefs)), nm[seq_along(final_coefs)], names(final_coefs))

# Keep only coefficients corresponding to PRS feature columns
coef_used <- final_coefs[intersect(names(final_coefs), colnames(all_prs_tune))]

Final_Coefficients <- data.frame(
  beta_matrix[, 1:4, drop = FALSE],
  BETA = as.matrix(beta_matrix[, names(coef_used), drop = FALSE]) %*% matrix(coef_used, ncol = 1)
)

# ---------------- Save coefficients ----------------
coef_file <- file.path(out_dir, sprintf("%s_Coefficients.csv", CV_number))
write.csv(Final_Coefficients, coef_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Coefficients saved to: %s\n", coef_file))

# ---------------- Predict ensemble rvPGS for tune/validation ----------------
predict_ens <- function(df_feats, coef_used) {
  as.numeric(as.matrix(df_feats[, names(coef_used), drop = TRUE]) %*% matrix(coef_used, ncol = 1))
}
PRS_Tune <- predict_ens(all_prs_tune, coef_used)
PRS_Val  <- predict_ens(all_prs_val,  coef_used)

PRS_Tune <- data.frame(IID = pheno_tune$IID, PRS = PRS_Tune)
PRS_Tune$scaled_PRS <- as.vector(scale(PRS_Tune$PRS))

PRS_Validation <- data.frame(IID = pheno_validation$IID, PRS = PRS_Val)
PRS_Validation$scaled_PRS <- as.vector(scale(PRS_Validation$PRS))

prs_tune_file <- file.path(out_dir, sprintf("%s_PRS_Tune.csv", CV_number))
prs_vali_file <- file.path(out_dir, sprintf("%s_PRS_Validation.csv", CV_number))
write.csv(PRS_Tune,       prs_tune_file, row.names = FALSE, quote = FALSE)
write.csv(PRS_Validation, prs_vali_file, row.names = FALSE, quote = FALSE)
cat(sprintf("Predictions written:\n  Tuning: %s\n  Validation: %s\n", prs_tune_file, prs_vali_file))

# ---------------- Diagnostics ----------------
cat("PRS_Tune scaled_PRS: range = ",
    paste(range(PRS_Tune$scaled_PRS, na.rm=TRUE), collapse=" - "),
    ", mean = ", mean(PRS_Tune$scaled_PRS, na.rm=TRUE),
    ", median = ", median(PRS_Tune$scaled_PRS, na.rm=TRUE), "\n")

cat("PRS_Validation scaled_PRS: range = ",
    paste(range(PRS_Validation$scaled_PRS, na.rm=TRUE), collapse=" - "),
    ", mean = ", mean(PRS_Validation$scaled_PRS, na.rm=TRUE),
    ", median = ", median(PRS_Validation$scaled_PRS, na.rm=TRUE), "\n")

# R^2 on residuals
pheno_tune$Ensemble_PRS       <- PRS_Tune$PRS
pheno_validation$Ensemble_PRS <- PRS_Validation$PRS

r2_tune <- summary(lm(y_tune ~ Ensemble_PRS,       data = pheno_tune))$r.squared
r2_val  <- summary(lm(y_validation ~ Ensemble_PRS, data = pheno_validation))$r.squared
cat(sprintf("R² (tuning residuals):     %.5f\n", r2_tune))
cat(sprintf("R² (validation residuals): %.5f\n", r2_val))

cat(sprintf("Successfully calculated ensemble rvPGS for CV %s at %s.\n",
            CV_number, Sys.time()))
