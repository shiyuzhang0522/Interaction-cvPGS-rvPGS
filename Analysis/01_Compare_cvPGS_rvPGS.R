########## Analysis to compare rvPGS and cvPGS ##########
########## Analysis/01_Compare_cvPGS_rvPGS.R ##########
########## Author: Shelley ##########
########## Date: 2025-08-08 ##########

# Working directory (kept as provided)
setwd("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Comparision")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(broom)
  library(nortest)    # ad.test
  library(tseries)    # jarque.bera.test
  library(moments)    # skewness, kurtosis
})

# ---------------- I/O ----------------
infile <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt"
data   <- fread(infile)

# Required columns (fail fast with a helpful message)
required_cols <- c("scaled_cvPRS","scaled_rvPRS","scaled_resid_TS_log_Zadj")
miss <- setdiff(required_cols, names(data))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# Coerce to numeric (in case anything read as character)
numify <- function(x) as.numeric(as.character(x))
data[, `:=`(
  scaled_cvPRS  = numify(scaled_cvPRS),
  scaled_rvPRS  = numify(scaled_rvPRS),
  scaled_resid_TS_log_Zadj = numify(scaled_resid_TS_log_Zadj)
)]

# ---------------- Module 1: Distribution (Hist + Density + Q-Q) ----------------
hist_density_plot <- function(df, col, color) {
  ggplot(df, aes(x = .data[[col]])) +
    geom_histogram(aes(y = after_stat(density)), fill = color, alpha = 0.4, bins = 50) +
    geom_density(color = color, linewidth = 1) +
    labs(x = paste("Scaled", col), y = "Density",
         title = paste("Density for", col, "(scaled)")) +
    theme_minimal(base_size = 12)
}

qq_plot <- function(df, col, color) {
  x <- df[[col]]
  x <- x[is.finite(x)]
  qq <- qqnorm(x, plot.it = FALSE)
  qq_df <- data.table(theoretical = qq$x, sample = qq$y)
  ggplot(qq_df, aes(x = theoretical, y = sample)) +
    geom_point(color = color, size = 1.2, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = color, linewidth = 1) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
         title = paste("Q–Q for", col, "(scaled)")) +
    theme_minimal(base_size = 12)
}

p_cv_hist <- hist_density_plot(data, "scaled_cvPRS", "#1f77b4")
p_cv_qq   <- qq_plot(data, "scaled_cvPRS", "#1f77b4")
p_rv_hist <- hist_density_plot(data, "scaled_rvPRS", "#d62728")
p_rv_qq   <- qq_plot(data, "scaled_rvPRS", "#d62728")

final_plot <- (p_cv_hist | p_cv_qq) / (p_rv_hist | p_rv_qq)
ggsave("cvPRS_rvPRS_distribution_QQ.pdf", final_plot, width = 10, height = 10)

# ---------------- Module 1b: Deviation from Normality ----------------
check_normality <- function(vec, name) {
  vec <- vec[is.finite(vec)]
  cat("\n---", name, "---\n")
  cat("n:", length(vec), "\n")
  cat("Skewness:", skewness(vec), "\n")
  cat("Kurtosis:", kurtosis(vec), "\n")
  # Tests (large n -> p ~ 0 even for tiny deviations; still useful for direction)
  cat("Kolmogorov–Smirnov p:",
      tryCatch(ks.test(vec, "pnorm", mean(vec), sd(vec))$p.value, error = function(e) NA), "\n")
  cat("Anderson–Darling p:",
      tryCatch(ad.test(vec)$p.value, error = function(e) NA), "\n")
  cat("Jarque–Bera p:",
      tryCatch(jarque.bera.test(vec)$p.value, error = function(e) NA), "\n")
}

check_normality(data$scaled_cvPRS, "scaled_cvPRS")
check_normality(data$scaled_rvPRS, "scaled_rvPRS")

# ---------------- Module 2: Effect sizes & Model fit ----------------
y  <- "scaled_resid_TS_log_Zadj"
X1 <- "scaled_cvPRS"
X2 <- "scaled_rvPRS"

# Fit models
model_cv   <- lm(as.formula(paste(y, "~", X1)), data = data)
model_rv   <- lm(as.formula(paste(y, "~", X2)), data = data)
model_both <- lm(as.formula(paste(y, "~", X1, "+", X2)), data = data)

tidy_with_r2 <- function(model, name) {
  tt <- as.data.table(broom::tidy(model))
  s  <- summary(model)
  tt[, `:=`(model = name, R2 = s$r.squared, adj_R2 = s$adj.r.squared)]
  tt[]
}

results_df <- rbindlist(list(
  tidy_with_r2(model_cv,   "cvPRS_only"),
  tidy_with_r2(model_rv,   "rvPRS_only"),
  tidy_with_r2(model_both, "both_PRS")
), fill = TRUE)

fwrite(results_df, "PRS_effect_size_and_model_fit.txt", sep = "\t")

# ΔR² using nested models
R2_cv   <- summary(model_cv)$r.squared
R2_rv   <- summary(model_rv)$r.squared
R2_both <- summary(model_both)$r.squared

adjR2_cv   <- summary(model_cv)$adj.r.squared
adjR2_rv   <- summary(model_rv)$adj.r.squared
adjR2_both <- summary(model_both)$adj.r.squared

anova_cv_then_rv <- anova(model_cv, model_both)  # add rv to cv
anova_rv_then_cv <- anova(model_rv, model_both)  # add cv to rv

delta_table <- data.table(
  Comparison = c("ΔR² (rvPRS | cvPRS)", "ΔR² (cvPRS | rvPRS)"),
  Delta_R2   = c(R2_both - R2_cv, R2_both - R2_rv),
  R2_base    = c(R2_cv, R2_rv),
  R2_full    = R2_both,
  AdjR2_base = c(adjR2_cv, adjR2_rv),
  AdjR2_full = adjR2_both,
  F_p_value  = c(anova_cv_then_rv$`Pr(>F)`[2], anova_rv_then_cv$`Pr(>F)`[2])
)
delta_table[, F_p_value := formatC(F_p_value, format = "e", digits = 2)]
fwrite(delta_table, "Delta_R2_summary.txt", sep = "\t", quote = FALSE)

# ---------------- Module 3: Mean ± SD across custom PRS percentiles ----------------
# Helper to compute stats across percentile bins (robust to ties)
compute_stats <- function(dt, prs_col, y_col,
                          breaks = c(0, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50,
                                     0.60, 0.70, 0.80, 0.90, 0.95, 0.99, 1),
                          labels = c("0-1%","1-5%","5-10%","10-20%","20-30%","30-40%",
                                     "40-50%","50-60%","60-70%","70-80%","80-90%",
                                     "90-95%","95-99%","99-100%")) {
  prs <- dt[[prs_col]]
  qs  <- quantile(prs, probs = breaks, na.rm = TRUE)
  qs  <- unique(qs) # collapse identical cut points if heavy ties
  b   <- cut(prs, breaks = qs, include.lowest = TRUE, ordered_result = TRUE)
  lvl <- levels(b)
  levels(b) <- if (length(lvl) == length(labels)) labels else labels[seq_along(lvl)]
  tmp <- data.table(bin = b, y = dt[[y_col]])
  out <- tmp[, .(mean_y = mean(y, na.rm = TRUE),
                 sd_y   = sd(y,   na.rm = TRUE),
                 n      = .N), by = bin]
  out[, bin_num := as.numeric(bin)]
  out[]
}

y_var <- "scaled_resid_TS_log_Zadj"
stats_cv <- compute_stats(data, "scaled_cvPRS", y_var)
stats_rv <- compute_stats(data, "scaled_rvPRS", y_var)

p_cv <- ggplot(stats_cv, aes(x = bin, y = mean_y, group = 1, color = bin_num)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.6) +
  geom_errorbar(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y), width = 0.15) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(x = "cvPRS percentile", y = y_var, title = "Mean ± SD of Y by cvPRS percentile") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p_rv <- ggplot(stats_rv, aes(x = bin, y = mean_y, group = 1, color = bin_num)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.6) +
  geom_errorbar(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y), width = 0.15) +
  scale_color_gradient(low = "mistyrose", high = "darkred") +
  labs(x = "rvPRS percentile", y = y_var, title = "Mean ± SD of Y by rvPRS percentile") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

final_plot2 <- p_cv / p_rv
ggsave("PRS_custom_percentile_mean_sd_vertical.pdf", final_plot2, width = 10, height = 8)
print(final_plot2)
