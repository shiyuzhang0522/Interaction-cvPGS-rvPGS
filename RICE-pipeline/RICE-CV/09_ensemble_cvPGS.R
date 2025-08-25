#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Rscript to get the ensemble cvPGS
# Author: Shelley
# Date: 2025-08-04
# Usage: Rscript 09_ensemble_cvPGS.R <CV_number>
# ------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages({
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
cat(sprintf("Processing for CV %s at %s.\n", CV_number, Sys.time()))

# ---------------- Base dir (edit this) ----------------
BASE_DIR <- "PATH/TO/PROJECT"

# Input dirs (keep filenames the same as your workflow)
CT_dir        <- file.path(BASE_DIR, "cvPGS/PRS_methods/C+T/C+T_PRS", CV_number)
LDpred2_dir   <- file.path(BASE_DIR, "LDpred2_LASSOsum2/LDpred2_LASSOsum2_estimates/PRS/LDpred2", CV_number)
LASSOsum2_dir <- file.path(BASE_DIR, "LDpred2_LASSOsum2/LDpred2_LASSOsum2_estimates/PRS/LASSOsum2", CV_number)

# Sample splits + phenotypes
sample_split_dir <- file.path(BASE_DIR, "Sample.Splitting")
pheno_file       <- file.path(sample_split_dir, "pheno_model.tsv")

# ---------------- Helpers ----------------
load_and_prefix <- function(file, prefix, key_id = c("#FID", "IID")) {
  dat <- fread(file, data.table = FALSE)
  if ("#FID" %in% names(dat)) colnames(dat)[colnames(dat) == "#FID"] <- "FID"
  non_id <- setdiff(colnames(dat), c("FID", "IID"))
  colnames(dat)[colnames(dat) %in% non_id] <- paste0(prefix, non_id)
  dat
}
load_rdata_and_prefix <- function(file, prefix, key_id = c("FID", "IID")) {
  e <- new.env()
  load(file, envir = e)
  obj <- e[[ls(e)[1]]]
  dat <- as.data.frame(obj)
  if ("#FID" %in% colnames(dat)) colnames(dat)[colnames(dat) == "#FID"] <- "FID"
  non_id <- setdiff(colnames(dat), c("FID", "IID"))
  colnames(dat)[colnames(dat) %in% non_id] <- paste0(prefix, non_id)
  dat
}

# ---------------- Read PRS matrices ----------------
# TRAIN
CT_prs_all_train     <- load_and_prefix(file.path(CT_dir, "C+T_prs_all_train.txt"), "CT_")
LDpred2_prs_train    <- load_rdata_and_prefix(file.path(LDpred2_dir, "LDpred2_train_prs_mat.RData"), "LDpred2_")
LASSOsum2_prs_train  <- load_rdata_and_prefix(file.path(LASSOsum2_dir, "LASSOsum2_train_prs_mat.RData"), "LASSOSum2_")

# TUNE
CT_prs_all_tune      <- load_and_prefix(file.path(CT_dir, "C+T_prs_all_tune.txt"), "CT_")
LDpred2_prs_tune     <- load_rdata_and_prefix(file.path(LDpred2_dir, "LDpred2_tune_prs_mat.RData"), "LDpred2_")
LASSOsum2_prs_tune   <- load_rdata_and_prefix(file.path(LASSOsum2_dir, "LASSOsum2_tune_prs_mat.RData"), "LASSOSum2_")

# VALIDATION
CT_prs_all_vali      <- load_and_prefix(file.path(CT_dir, "C+T_prs_all_validation.txt"), "CT_")
LDpred2_prs_vali     <- load_rdata_and_prefix(file.path(LDpred2_dir, "LDpred2_vali_prs_mat.RData"), "LDpred2_")
LASSOsum2_prs_vali   <- load_rdata_and_prefix(file.path(LASSOsum2_dir, "LASSOsum2_vali_prs_mat.RData"), "LASSOSum2_")

# Merge by FID/IID
prs_train_all <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID")),
                        list(CT_prs_all_train, LDpred2_prs_train, LASSOsum2_prs_train))
prs_tune_all <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID")),
                       list(CT_prs_all_tune, LDpred2_prs_tune, LASSOsum2_prs_tune))
prs_validation_all <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID")),
                             list(CT_prs_all_vali, LDpred2_prs_vali, LASSOsum2_prs_vali))
cat(sprintf("Read all PRS matrices for CV %s at %s.\n", CV_number, Sys.time()))

# ---------------- Phenotypes & splits ----------------
pheno_all <- fread(pheno_file, data.table = FALSE)

train_ids <- fread(file.path(sample_split_dir, sprintf("cv%s_train.txt", CV_number)), header = FALSE)[[1]]
tune_ids  <- fread(file.path(sample_split_dir, sprintf("cv%s_tune.txt",  CV_number)), header = FALSE)[[1]]
vali_ids  <- fread(file.path(sample_split_dir, sprintf("cv%s_val.txt",   CV_number)), header = FALSE)[[1]]

pheno_train      <- pheno_all %>% dplyr::filter(eid %in% train_ids) %>% dplyr::mutate(IID = eid)
pheno_tune       <- pheno_all %>% dplyr::filter(eid %in% tune_ids)  %>% dplyr::mutate(IID = eid)
pheno_validation <- pheno_all %>% dplyr::filter(eid %in% vali_ids)  %>% dplyr::mutate(IID = eid)

pheno_train      <- dplyr::left_join(pheno_train,      prs_train_all,      by = c("IID" = "IID", "eid" = "FID"))
pheno_tune       <- dplyr::left_join(pheno_tune,       prs_tune_all,       by = c("IID" = "IID", "eid" = "FID"))
pheno_validation <- dplyr::left_join(pheno_validation, prs_validation_all, by = c("IID" = "IID", "eid" = "FID"))

# ---------------- Select PRS columns ----------------
prs_col_prefix <- c("CT_", "LDpred2_", "LASSOSum2_")
prs_cols <- grep(paste0("^(", paste(prs_col_prefix, collapse = "|"), ")"), colnames(pheno_train), value = TRUE)

prs_train_all       <- pheno_train[, prs_cols, drop = FALSE]
prs_tune_all        <- pheno_tune[, prs_cols, drop = FALSE]
prs_validation_all  <- pheno_validation[, prs_cols, drop = FALSE]

cat(sprintf("PRS Train: %d rows, %d columns\n", nrow(prs_train_all), ncol(prs_train_all)))
cat(sprintf("PRS Tune: %d rows, %d columns\n", nrow(prs_tune_all), ncol(prs_tune_all)))
cat(sprintf("PRS Validation: %d rows, %d columns\n", nrow(prs_validation_all), ncol(prs_validation_all)))

# ---------------- Drop collinear / highly correlated features ----------------
drop <- caret::findLinearCombos(prs_tune_all)$remove
drop <- names(prs_tune_all)[drop]
cat("Perfect linear combinations to drop:\n"); print(drop)

prs_train_all      <- dplyr::select(prs_train_all,      -dplyr::all_of(drop))
prs_tune_all       <- dplyr::select(prs_tune_all,       -dplyr::all_of(drop))
prs_validation_all <- dplyr::select(prs_validation_all, -dplyr::all_of(drop))

mtx  <- stats::cor(prs_tune_all)
drop <- caret::findCorrelation(mtx, cutoff = 0.98)
drop <- names(prs_tune_all)[drop]
cat("Highly correlated (>0.98) to drop:\n"); print(drop)

prs_train_all      <- dplyr::select(prs_train_all,      -dplyr::all_of(drop))
prs_tune_all       <- dplyr::select(prs_tune_all,       -dplyr::all_of(drop))
prs_validation_all <- dplyr::select(prs_validation_all, -dplyr::all_of(drop))

cat("Dimensions after filtering:\n")
cat("prs_train_all: ", dim(prs_train_all)[1], "x", dim(prs_train_all)[2], "\n")
cat("prs_tune_all: ", dim(prs_tune_all)[1],  "x", dim(prs_tune_all)[2],  "\n")
cat("prs_validation_all: ", dim(prs_validation_all)[1], "x", dim(prs_validation_all)[2], "\n")

# ---------------- Ensemble functions ----------------
Ensemble_Function_Continuous <- function(x,y){
  x <- as.matrix(x[!is.na(y),])
  y <- y[!is.na(y)]
  lasso_mod <- glmnet::cv.glmnet(x, y, family = "gaussian", alpha = 1, type.measure = "mse", nfolds = 10)
  ridge_mod <- glmnet::cv.glmnet(x, y, family = "gaussian", alpha = 0, type.measure = "mse", nfolds = 10)
  lasso_prediction_x <- predict(lasso_mod, x)
  ridge_prediction_x <- predict(ridge_mod, x)
  ensemble_mod <- lm(y ~ ., data = data.frame(lasso_prediction_x, ridge_prediction_x))
  ensemble_prediction_x <- ensemble_mod$fitted
  coefficients_x <- coef(lm(y ~ ., data.frame(y = ensemble_prediction_x, x)))
  list(Coefficients = coefficients_x)
}
Ensemble_Function_Binary <- function(x,y){
  x <- as.matrix(x[!is.na(y),])
  y <- y[!is.na(y)]
  lasso_mod <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "auc")
  ridge_mod <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 0, type.measure = "auc")
  lasso_prediction_x <- predict(lasso_mod, x, type = "link")
  ridge_prediction_x <- predict(ridge_mod, x, type = "link")
  ensemble_mod <- glm(y ~ ., data = data.frame(lasso_prediction_x, ridge_prediction_x), family = binomial())
  ensemble_prediction_x <- predict(ensemble_mod, data.frame(lasso_prediction_x, ridge_prediction_x), type = "link")
  coefficients_x <- coef(lm(y ~ ., data.frame(y = ensemble_prediction_x, x)))
  list(Coefficients = coefficients_x)
}
Ensemble_Function <- function(x,y,family = c("continuous","binary")){
  family <- match.arg(family)
  if (family == "continuous") Ensemble_Function_Continuous(x,y) else Ensemble_Function_Binary(x,y)
}

# ---------------- Residualize TL (null models) ----------------
model.null <- lm(TS_log_Zadj ~ age + age2 + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_train)
y_train <- model.null$residuals

model.null <- lm(TS_log_Zadj ~ age + age2 + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_tune)
y_tune <- model.null$residuals

model.null <- lm(TS_log_Zadj ~ age + age2 + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = pheno_validation)
y_validation <- model.null$residuals

# Sanity prints
cat("Length y_train / n(pheno_train):", length(y_train), nrow(pheno_train), "\n")
cat("Length y_tune  / n(pheno_tune):",  length(y_tune),  nrow(pheno_tune),  "\n")
cat("Length y_valid / n(pheno_valid):", length(y_validation), nrow(pheno_validation), "\n")

# Attach residuals
pheno_tune$y_tune <- NA;        pheno_tune$y_tune[!is.na(pheno_tune$TS_log_Zadj)] <- y_tune
pheno_validation$y_validation <- NA; pheno_validation$y_validation[!is.na(pheno_validation$TS_log_Zadj)] <- y_validation

# ---------------- Fit ensemble on tuning ----------------
Results <- Ensemble_Function(x = prs_tune_all, y = pheno_tune[,"y_tune"], family = "continuous")
Results$Coefficients[is.na(Results$Coefficients)] <- 0

# ---------------- Save coefficients ----------------
out_coef_dir <- file.path(BASE_DIR, "cvPGS/PRS_methods/Ensemble.cvPGS", CV_number)
dir.create(out_coef_dir, recursive = TRUE, showWarnings = FALSE)
coefs <- data.frame(feature = rownames(Results$Coefficients), Coefficient = as.numeric(Results$Coefficients))
write.csv(coefs, file.path(out_coef_dir, sprintf("CV%s_Coefficients.csv", CV_number)), row.names = FALSE, quote = FALSE)
cat(sprintf("Coefficients saved to: %s\n", file.path(out_coef_dir, sprintf("CV%s_Coefficients.csv", CV_number))))

# ---------------- Compute ensemble PRS for train/tune/validation ----------------
feat <- names(Results$Coefficients)[-1]
PRS_Train      <- as.matrix(pheno_train[,      feat]) %*% matrix(Results$Coefficients[-1], ncol = 1)
PRS_Tune       <- as.matrix(pheno_tune[,       feat]) %*% matrix(Results$Coefficients[-1], ncol = 1)
PRS_Validation <- as.matrix(pheno_validation[, feat]) %*% matrix(Results$Coefficients[-1], ncol = 1)

out_dir <- file.path(BASE_DIR, "cvPGS/PRS_methods/Ensemble.cvPGS", CV_number)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_prs <- function(scores, ids, path) {
  df <- data.frame(IID = ids, PRS = as.vector(scores))
  df$scaled_PRS <- as.vector(scale(df$PRS))
  write.csv(df, path, row.names = FALSE, quote = FALSE)
}
save_prs(PRS_Train,      pheno_train$IID,      file.path(out_dir, "Ensemble_PRS_Train.csv"))
save_prs(PRS_Tune,       pheno_tune$IID,       file.path(out_dir, "Ensemble_PRS_Tune.csv"))
save_prs(PRS_Validation, pheno_validation$IID, file.path(out_dir, "Ensemble_PRS_Validation.csv"))

cat("PRS_Train scaled_PRS: range = ", paste(range(scale(PRS_Train), na.rm=TRUE), collapse=" - "),
    ", mean = ", mean(scale(PRS_Train), na.rm=TRUE),
    ", median = ", median(scale(PRS_Train), na.rm=TRUE), "\n")
cat("PRS_Tune  scaled_PRS: range = ", paste(range(scale(PRS_Tune), na.rm=TRUE), collapse=" - "),
    ", mean = ", mean(scale(PRS_Tune), na.rm=TRUE),
    ", median = ", median(scale(PRS_Tune), na.rm=TRUE), "\n")
cat("PRS_Valid scaled_PRS: range = ", paste(range(scale(PRS_Validation), na.rm=TRUE), collapse=" - "),
    ", mean = ", mean(scale(PRS_Validation), na.rm=TRUE),
    ", median = ", median(scale(PRS_Validation), na.rm=TRUE), "\n")

# ---------------- R² on residuals ----------------
pheno_validation$Ensemble_PRS <- as.vector(PRS_Validation)
r2_val  <- summary(lm(y_validation ~ Ensemble_PRS, data = pheno_validation))$r.squared
pheno_tune$Ensemble_PRS <- as.vector(PRS_Tune)
r2_tune <- summary(lm(y_tune ~ Ensemble_PRS, data = pheno_tune))$r.squared
cat(sprintf("R² (tune residuals): %.5f\n", r2_tune))
cat(sprintf("R² (validation residuals): %.5f\n", r2_val))

cat(sprintf("Successfully generated the ensemble cvPGS for CV %s at %s.\n", CV_number, Sys.time()))
cat("Ensemble PRS coefficients and predictions saved successfully.\n")
