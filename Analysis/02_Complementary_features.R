#### 02. Explore the complementary features of cvPGS and rvPGS ##################
############## Rscript to explore the complementary features of cvPGS and rvPGS ##############
############## Author: Shelley #############
############## Date: 2025-08-08 #############

# --- Setup --------------------------------------------------------------------
setwd("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Comparision")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(broom)
})

# --- Input --------------------------------------------------------------------
data <- fread("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt")

# Require columns
req <- c("scaled_cvPRS", "scaled_rvPRS", "scaled_resid_TS_log_Zadj")
miss <- setdiff(req, names(data))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# Ensure numeric
numify <- function(x) as.numeric(as.character(x))
data[, `:=`(
  scaled_cvPRS = numify(scaled_cvPRS),
  scaled_rvPRS = numify(scaled_rvPRS),
  scaled_resid_TS_log_Zadj = numify(scaled_resid_TS_log_Zadj)
)]

dt <- copy(data)

# --- Phenotype highlight groups ------------------------------------------------
# 1) TL_decile_group: top/bottom 10% of phenotype
ph_top10 <- quantile(dt$scaled_resid_TS_log_Zadj, 0.90, na.rm = TRUE)
ph_bot10 <- quantile(dt$scaled_resid_TS_log_Zadj, 0.10, na.rm = TRUE)
dt[, TL_decile_group := fifelse(
  scaled_resid_TS_log_Zadj >= ph_top10, "Highest 10%",
  fifelse(scaled_resid_TS_log_Zadj <= ph_bot10, "Lowest 10%", "Other")
)]

# 2) TL_outlier_group: ±3 SD on phenotype
dt[, TL_outlier_group := fifelse(
  scaled_resid_TS_log_Zadj >= 3,  "Highest outliers",
  fifelse(scaled_resid_TS_log_Zadj <= -3, "Lowest outliers", "other")
)]

# --- PRS thresholds (centralized definitions) ---------------------------------
prs_thresholds <- list(
  sd3 = list(
    cv_low = function(z) -3, cv_hi = function(z) 3,
    rv_low = function(z) -3, rv_hi = function(z) 3,
    label  = "±3 SD"
  ),
  top1 = list(
    cv_low = function(z) quantile(z, 0.01, na.rm = TRUE),
    cv_hi  = function(z) quantile(z, 0.99, na.rm = TRUE),
    rv_low = function(z) quantile(z, 0.01, na.rm = TRUE),
    rv_hi  = function(z) quantile(z, 0.99, na.rm = TRUE),
    label  = "top/bottom 1%"
  ),
  top5 = list(
    cv_low = function(z) quantile(z, 0.05, na.rm = TRUE),
    cv_hi  = function(z) quantile(z, 0.95, na.rm = TRUE),
    rv_low = function(z) quantile(z, 0.05, na.rm = TRUE),
    rv_hi  = function(z) quantile(z, 0.95, na.rm = TRUE),
    label  = "top/bottom 5%"
  ),
  top10 = list(
    cv_low = function(z) quantile(z, 0.10, na.rm = TRUE),
    cv_hi  = function(z) quantile(z, 0.90, na.rm = TRUE),
    rv_low = function(z) quantile(z, 0.10, na.rm = TRUE),
    rv_hi  = function(z) quantile(z, 0.90, na.rm = TRUE),
    label  = "top/bottom 10%"
  )
)
threshold_order <- c("sd3","top1","top5","top10")

# --- Helper: summarize identification under a threshold ------------------------
summarize_identification <- function(df, pheno_group_col, th_spec) {
  df <- copy(df)

  # Cutoffs
  cv_low <- th_spec$cv_low(df$scaled_cvPRS);  cv_hi <- th_spec$cv_hi(df$scaled_cvPRS)
  rv_low <- th_spec$rv_low(df$scaled_rvPRS);  rv_hi <- th_spec$rv_hi(df$scaled_rvPRS)

  # PRS flags
  df[, cv_high   := scaled_cvPRS >= cv_hi]
  df[, cv_low_gp := scaled_cvPRS <= cv_low]
  df[, rv_high   := scaled_rvPRS >= rv_hi]
  df[, rv_low_gp := scaled_rvPRS <= rv_low]

  # Phenotype totals
  total_high <- df[get(pheno_group_col) %like% "Highest", .N]
  total_low  <- df[get(pheno_group_col) %like% "Lowest",  .N]

  # Captures
  cv_high_id <- df[get(pheno_group_col) %like% "Highest" & cv_high   == TRUE, .N]
  cv_low_id  <- df[get(pheno_group_col) %like% "Lowest"  & cv_low_gp == TRUE, .N]
  rv_high_id <- df[get(pheno_group_col) %like% "Highest" & rv_high   == TRUE, .N]
  rv_low_id  <- df[get(pheno_group_col) %like% "Lowest"  & rv_low_gp == TRUE, .N]

  # Sizes of PRS groups
  cv_high_N <- df[cv_high   == TRUE, .N]
  cv_low_N  <- df[cv_low_gp == TRUE, .N]
  rv_high_N <- df[rv_high   == TRUE, .N]
  rv_low_N  <- df[rv_low_gp == TRUE, .N]

  safe_prop <- function(num, den) ifelse(den > 0, num / den, NA_real_)

  data.table(
    Threshold = th_spec$label,
    Group     = c("High", "Low"),
    Total     = c(total_high, total_low),

    cvPRS_identified_N   = c(cv_high_id, cv_low_id),
    cvPRS_identified_pct = round(c(cv_high_id/total_high, cv_low_id/total_low) * 100, 2),
    rvPRS_identified_N   = c(rv_high_id, rv_low_id),
    rvPRS_identified_pct = round(c(rv_high_id/total_high, rv_low_id/total_low) * 100, 2),

    cvPRS_threshold_N = c(cv_high_N, cv_low_N),
    rvPRS_threshold_N = c(rv_high_N, rv_low_N),

    cvPRS_identified_prop_of_threshold = round(
      c(safe_prop(cv_high_id, cv_high_N), safe_prop(cv_low_id, cv_low_N)) * 100, 2),
    rvPRS_identified_prop_of_threshold = round(
      c(safe_prop(rv_high_id, rv_high_N), safe_prop(rv_low_id, rv_low_N)) * 100, 2)
  )
}

# --- Identification summaries --------------------------------------------------
summary_decile <- rbindlist(lapply(prs_thresholds, function(th) {
  summarize_identification(copy(dt), "TL_decile_group", th)
}), use.names = TRUE, fill = TRUE)

summary_outlier <- rbindlist(lapply(prs_thresholds, function(th) {
  summarize_identification(copy(dt), "TL_outlier_group", th)
}), use.names = TRUE, fill = TRUE)

# Save
fwrite(summary_decile,  "complementary.identification.summary_phenotype_decile.txt",           sep = "\t", quote = FALSE)
fwrite(summary_outlier, "complementary.identification.summary_phenotype_summary_outlier.txt", sep = "\t", quote = FALSE)

# --- Scatter plots: who identifies TL extremes (±3 SD PRS) --------------------
build_sd3_plot_byPRS <- function(d, tl_col,
                                 title_prefix = "",
                                 col_cv = "#0571b0",    # blue
                                 col_rv = "#ca0020",    # red
                                 col_both = "#6a3d9a"   # purple
) {
  dd <- copy(d)

  # ±3 SD PRS flags
  dd[, cv_hi := scaled_cvPRS >=  3]
  dd[, cv_lo := scaled_cvPRS <= -3]
  dd[, rv_hi := scaled_rvPRS >=  3]
  dd[, rv_lo := scaled_rvPRS <= -3]

  # Phenotype high/low flags
  dd[, TL_high := (get(tl_col) %like% "Highest")]
  dd[, TL_low  := (get(tl_col) %like% "Lowest")]

  # Identified-by flags
  dd[, TL_high_identified := TL_high & (cv_hi | rv_hi)]
  dd[, TL_low_identified  := TL_low  & (cv_lo | rv_lo)]

  dd[, src_high := fifelse(TL_high_identified & cv_hi &  rv_hi, "both",
                    fifelse(TL_high_identified & cv_hi & !rv_hi, "cv_only",
                    fifelse(TL_high_identified & !cv_hi & rv_hi, "rv_only", NA_character_)))]
  dd[, src_low  := fifelse(TL_low_identified  & cv_lo &  rv_lo, "both",
                    fifelse(TL_low_identified  & cv_lo & !rv_lo, "cv_only",
                    fifelse(TL_low_identified  & !cv_lo & rv_lo, "rv_only", NA_character_)))]

  dd[, identified_by := NA_character_]
  dd[src_high == "cv_only" | src_low == "cv_only", identified_by := "cvPGS"]
  dd[src_high == "rv_only" | src_low == "rv_only", identified_by := "rvPGS"]
  dd[src_high == "both"    | src_low == "both",    identified_by := "both"]
  dd[, identified_by := factor(identified_by, levels = c("cvPGS","rvPGS","both"))]

  # Counts for annotation
  cnt <- list(
    high_cv_only = sum(dd$src_high == "cv_only", na.rm = TRUE),
    high_rv_only = sum(dd$src_high == "rv_only", na.rm = TRUE),
    high_both    = sum(dd$src_high == "both",    na.rm = TRUE),
    high_total   = sum(dd$TL_high_identified,    na.rm = TRUE),
    low_cv_only  = sum(dd$src_low  == "cv_only", na.rm = TRUE),
    low_rv_only  = sum(dd$src_low  == "rv_only", na.rm = TRUE),
    low_both     = sum(dd$src_low  == "both",    na.rm = TRUE),
    low_total    = sum(dd$TL_low_identified,     na.rm = TRUE)
  )
  annot_label <- sprintf(
    "High: cv-only=%d  rv-only=%d  both=%d  total=%d\nLow:  cv-only=%d  rv-only=%d  both=%d  total=%d",
    cnt$high_cv_only, cnt$high_rv_only, cnt$high_both, cnt$high_total,
    cnt$low_cv_only,  cnt$low_rv_only,  cnt$low_both,  cnt$low_total
  )

  ggplot(dd, aes(scaled_cvPRS, scaled_rvPRS)) +
    geom_point(data = dd[is.na(identified_by)], color = "grey80", alpha = 0.4, size = 0.9) +
    geom_point(data = dd[!is.na(identified_by)], aes(color = identified_by), size = 1.6, alpha = 0.9) +
    geom_vline(xintercept = c(-3, 3), linetype = "dashed") +
    geom_hline(yintercept = c(-3, 3), linetype = "dashed") +
    scale_color_manual(values = c("cvPGS" = col_cv, "rvPGS" = col_rv, "both" = col_both),
                       limits = levels(dd$identified_by), drop = FALSE) +
    labs(x = "cvPGS (scaled)", y = "rvPGS (scaled)", color = "Identified by",
         title = paste0(title_prefix, "\nPRS thresholds: ±3 SD")) +
    annotate("text", x = Inf, y = Inf, hjust = 1.02, vjust = 1.3, size = 3.2, label = annot_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

p_decile  <- build_sd3_plot_byPRS(dt, "TL_decile_group",
                                  "TL extremes (top/bottom 10% of phenotype)")
p_outlier <- build_sd3_plot_byPRS(dt, "TL_outlier_group",
                                  "TL outliers (|scaled TL| ≥ 3)")

# --- Enrichment: phenotype outliers within PRS outliers ------------------------
# Combined phenotype outlier flag (either highest or lowest)
dt[, pheno_outlier := TL_outlier_group %in% c("Highest outliers","Lowest outliers")]

calc_enrichment_combined <- function(df, th_key, th_spec) {
  d <- copy(df)
  # PRS cutoffs
  cv_lo <- th_spec$cv_low(d$scaled_cvPRS);  cv_hi <- th_spec$cv_hi(d$scaled_cvPRS)
  rv_lo <- th_spec$rv_low(d$scaled_rvPRS);  rv_hi <- th_spec$rv_hi(d$scaled_rvPRS)
  # PRS outlier flags (either tail)
  d[, cv_out := (scaled_cvPRS <= cv_lo | scaled_cvPRS >= cv_hi)]
  d[, rv_out := (scaled_rvPRS <= rv_lo | scaled_rvPRS >= rv_hi)]

  fe_row <- function(out_flag, prs_name) {
    tab <- table(PRS_outlier = out_flag, Pheno_outlier = d$pheno_outlier)
    if (!all(dim(tab) == c(2,2))) {
      full <- matrix(0, nrow=2, ncol=2,
                     dimnames=list(PRS_outlier=c(FALSE,TRUE), Pheno_outlier=c(FALSE,TRUE)))
      full[rownames(tab), colnames(tab)] <- tab
      tab <- full
    }
    ft <- fisher.test(tab)
    data.table(
      Threshold_key = th_key,
      Threshold     = th_spec$label,
      PRS_type      = prs_name,
      N_total       = nrow(d),
      N_PRS_group   = sum(tab[2, ]),
      N_pheno_out   = sum(tab[, 2]),
      N_both        = tab[2, 2],
      OR            = unname(ft$estimate),
      CI_low        = ft$conf.int[1],
      CI_high       = ft$conf.int[2],
      p_value       = ft$p.value
    )
  }

  rbindlist(list(
    fe_row(d$cv_out, "cvPGS"),
    fe_row(d$rv_out, "rvPGS")
  ))
}

enrich_combined <- rbindlist(lapply(threshold_order, function(k) {
  calc_enrichment_combined(dt, k, prs_thresholds[[k]])
}), use.names = TRUE)

# Pretty factors and save
enrich_combined[, Threshold := factor(Threshold,
  levels = c("±3 SD","top/bottom 1%","top/bottom 5%","top/bottom 10%"))]
enrich_combined[, PRS_type := factor(PRS_type, levels = c("cvPGS","rvPGS"))]
fwrite(enrich_combined, "Enrichment_TLoutliers_combined.tsv", sep = "\t", quote = FALSE)

# Plot: OR + 95% CI (log scale)
p_combined <- ggplot(enrich_combined,
                     aes(x = Threshold, y = OR, ymin = CI_low, ymax = CI_high, color = PRS_type)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.5), linewidth = 0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("cvPGS" = "#1f77b4", "rvPGS" = "#d62728"), name = "PRS") +
  labs(
    title = "Enrichment of phenotype outliers within PRS outlier groups (combined high + low)",
    subtitle = "Odds ratio with 95% CI (log scale); PRS outlier = either tail (≤ low OR ≥ high)",
    x = "PRS threshold", y = "Odds ratio (log)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# --- Final figure: identification panels + enrichment --------------------------
top_row <- p_decile | p_outlier
final_figure <- (top_row / p_combined) +
  plot_layout(heights = c(1, 1.1), guides = "collect") &
  theme(legend.position = "right")

ggsave("PRS_Complementary_plots.pdf", final_figure, width = 12, height = 12)
print(final_figure)

# --- Also save identification summaries used earlier ---------------------------
# (kept for compatibility with your original script)
fwrite(summary_decile,  "complementary.identification.summary_phenotype_decile.txt",           sep = "\t", quote = FALSE)
fwrite(summary_outlier, "complementary.identification.summary_phenotype_summary_outlier.txt", sep = "\t", quote = FALSE)
