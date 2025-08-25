#### 03. Correlation between scaled_cvPRS and scaled_rvPRS ######################
######## Rscript to explore the correlation between scaled_cvPRS and scaled_rvPRS
######## Author: Shelley ########
######## Date: 2025-08-08 ########

# --- Setup --------------------------------------------------------------------
setwd("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Correlations")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(broom)
})

# --- Input --------------------------------------------------------------------
data <- fread("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt")

# --- 1) Linear model ----------------------------------------------------------
lm_fit <- lm(scaled_rvPRS ~ scaled_cvPRS, data = data)
summary_fit <- summary(lm_fit)

lm_results <- tidy(lm_fit)  # broom output
print(lm_results)

# Save
fwrite(lm_results, "lm_results_cvPRS_vs_rvPRS.txt", sep = "\t")

# --- 2) Pearson correlation ---------------------------------------------------
cor_test <- cor.test(data$scaled_cvPRS, data$scaled_rvPRS, use = "complete.obs")
print(cor_test)

# --- 3) Stratified analysis by cvPRS tails ------------------------------------
dt <- as.data.table(data)

analyze_subset <- function(D, label) {
  D <- D[complete.cases(scaled_cvPRS, scaled_rvPRS)]
  if (nrow(D) < 3) {
    return(data.table(
      Subset = label, N = nrow(D),
      Pearson_r = NA_real_, r_CI_low = NA_real_, r_CI_high = NA_real_, r_p = NA_real_,
      LM_slope = NA_real_, slope_CI_low = NA_real_, slope_CI_high = NA_real_, slope_p = NA_real_,
      R2 = NA_real_
    ))
  }
  # Pearson
  ct <- cor.test(D$scaled_cvPRS, D$scaled_rvPRS, method = "pearson", conf.level = 0.95)
  # Linear model
  fit <- lm(scaled_rvPRS ~ scaled_cvPRS, data = D)
  sm  <- summary(fit)
  ci  <- confint(fit, level = 0.95)

  data.table(
    Subset       = label,
    N            = nrow(D),
    Pearson_r    = unname(ct$estimate),
    r_CI_low     = unname(ct$conf.int[1]),
    r_CI_high    = unname(ct$conf.int[2]),
    r_p          = ct$p.value,
    LM_slope     = unname(coef(fit)["scaled_cvPRS"]),
    slope_CI_low = unname(ci["scaled_cvPRS", 1]),
    slope_CI_high= unname(ci["scaled_cvPRS", 2]),
    slope_p      = unname(coef(sm)["scaled_cvPRS", "Pr(>|t|)"]),
    R2           = sm$r.squared
  )
}

# Define thresholds
q10 <- quantile(dt$scaled_cvPRS, 0.90, na.rm = TRUE)
q05 <- quantile(dt$scaled_cvPRS, 0.95, na.rm = TRUE)
q01 <- quantile(dt$scaled_cvPRS, 0.99, na.rm = TRUE)

res_list <- list(
  analyze_subset(dt,                      "All"),
  analyze_subset(dt[scaled_cvPRS >= q10], "Top 10% (cvPRS)"),
  analyze_subset(dt[scaled_cvPRS >= q05], "Top 5% (cvPRS)"),
  analyze_subset(dt[scaled_cvPRS >= q01], "Top 1% (cvPRS)")
)

cvprs_toptail_corr <- rbindlist(res_list, use.names = TRUE)

# Save
fwrite(cvprs_toptail_corr, "cvPRS_rvPRS_correlation_topTails_withCIs.txt", sep = "\t")
print(cvprs_toptail_corr)

# --- 4) Plot Pearson r with 95% CI --------------------------------------------
cvprs_toptail_corr[, Subset := factor(
  Subset,
  levels = c("All", "Top 10% (cvPRS)", "Top 5% (cvPRS)", "Top 1% (cvPRS)")
)]

p_r_ci <- ggplot(cvprs_toptail_corr,
                 aes(x = Subset, y = Pearson_r, ymin = r_CI_low, ymax = r_CI_high)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(linewidth = 1) +
  geom_point(size = 2) +
  geom_text(aes(label = paste0("N=", N), y = r_CI_high),
            vjust = -0.6, size = 3) +
  labs(
    title = "Correlation (Pearson r) between cvPRS and rvPRS by cvPRS tail",
    x = NULL, y = "Pearson r (95% CI)"
  ) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_r_ci)
ggsave("cvPRS_rvPRS_Pearson_r_CIs_by_tail.pdf", p_r_ci, width = 4, height = 4)
