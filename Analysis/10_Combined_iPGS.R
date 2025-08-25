################ Rscript 10: Combine cvPGS + rvPGS & assess prediction ################
################ Author: Shelley #######################################################
################ Date: 2025-08-19 #####################################################

# Working directory
setwd("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Combined.RICE.PRS")

# ---- Packages -----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(broom)
  library(ggpubr)
  library(pROC)
})

# ---- Data ---------------------------------------------------------------------------
infile <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt"
dt <- fread(infile)

# Keep only columns we need and complete cases
set.seed(20250819)
y_col <- "scaled_resid_TS_log_Zadj"  # outcome already residualized
cv_col <- "scaled_cvPRS"
rv_col <- "scaled_rvPRS"

need <- c(y_col, cv_col, rv_col)
d <- dt[complete.cases(dt[, ..need]), ..need]

# ---- Train/test split (70/30) -------------------------------------------------------
idx <- sample(seq_len(nrow(d)), size = floor(0.7 * nrow(d)))
train <- d[idx]
test  <- d[-idx]

# ---- Models (linear) ----------------------------------------------------------------
f_cv  <- as.formula(paste(y_col, "~", cv_col))
f_rv  <- as.formula(paste(y_col, "~", rv_col))
f_add <- as.formula(paste(y_col, "~", cv_col, "+", rv_col))
f_int <- as.formula(paste(y_col, "~", cv_col, "*", rv_col))

m_cv  <- lm(f_cv,  data = train)
m_rv  <- lm(f_rv,  data = train)
m_add <- lm(f_add, data = train)
m_int <- lm(f_int, data = train)

# ---- Test-set metrics ---------------------------------------------------------------
metrics <- function(obs, pred) {
  sse <- sum((obs - pred)^2)
  sst <- sum((obs - mean(obs))^2)
  data.frame(
    R2   = 1 - sse/sst,
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE  = mean(abs(obs - pred))
  )
}

# Predictions
obs      <- test[[y_col]]
pred_cv  <- predict(m_cv,  newdata = test)
pred_rv  <- predict(m_rv,  newdata = test)
pred_add <- predict(m_add, newdata = test)
pred_int <- predict(m_int, newdata = test)

# Performance table
perf <- bind_rows(
  cbind(model = "cvPRS only",      metrics(obs, pred_cv)),
  cbind(model = "rvPRS only",      metrics(obs, pred_rv)),
  cbind(model = "additive.PRS",    metrics(obs, pred_add)),
  cbind(model = "interaction.PRS", metrics(obs, pred_int))
) %>% as.data.frame()

# Deltas vs additive
perf$delta_R2_vs_additive <- perf$R2 - perf$R2[perf$model == "additive.PRS"]

# P-values vs additive via paired t-test on squared errors
se_add <- (obs - pred_add)^2
p_vs_add <- function(pred_other) {
  se_other <- (obs - pred_other)^2
  t.test(se_other, se_add, paired = TRUE)$p.value
}
perf$p_vs_additive <- NA_real_
perf$p_vs_additive[perf$model == "cvPRS only"]      <- p_vs_add(pred_cv)
perf$p_vs_additive[perf$model == "rvPRS only"]      <- p_vs_add(pred_rv)
perf$p_vs_additive[perf$model == "interaction.PRS"] <- p_vs_add(pred_int)

# ---- Save performance table + a publication-style table figure ----------------------
fwrite(perf, file = "PRS_test_performance_metrics.txt", sep = "\t")

# Pretty table (as PDF)
order_models <- c("cvPRS only","rvPRS only","additive.PRS","interaction.PRS")
perf_plot <- perf %>%
  mutate(model = factor(model, levels = order_models)) %>%
  arrange(model) %>%
  mutate(across(where(is.numeric), ~ format(.x, scientific = TRUE, digits = 5)))

my_theme <- ttheme(
  base_style = "classic",
  colnames.style = colnames_style(fill = "grey90", color = "black", face = "bold", size = 12),
  tbody.style    = tbody_style(color = "black", size = 11, hjust = 0.5, fill = c("white", "grey98"))
)
tbl <- ggtexttable(perf_plot, rows = NULL, theme = my_theme)
ggexport(tbl, filename = "PRS_test_performance_table.pdf", width = 9, height = 3)

# ---- Extremes classification: AUCs for cvPRS/rvPRS/interaction ----------------------
# Helper labelers
mk_labels_pct <- function(y, p, side = c("top","bottom")) {
  side <- match.arg(side)
  cut <- if (side == "top") quantile(y, 1 - p, na.rm = TRUE) else quantile(y, p, na.rm = TRUE)
  as.integer(if (side == "top") y >= cut else y <= cut)
}
mk_labels_out <- function(y, kind = c("high","low")) {
  kind <- match.arg(kind)
  as.integer(if (kind == "high") y > 3 else y < -3)
}
auc_ci <- function(labels, score) {
  if (length(unique(labels)) < 2L) return(c(NA, NA, NA))
  r <- pROC::roc(labels, score, quiet = TRUE, direction = "<")
  ci <- suppressMessages(pROC::ci.auc(r, conf.level = 0.95, method = "delong"))
  c(as.numeric(r$auc), as.numeric(ci[1]), as.numeric(ci[3]))
}
delong_p <- function(labels, s1, s2) {
  if (length(unique(labels)) < 2L) return(NA_real_)
  r1 <- tryCatch(pROC::roc(labels, s1, quiet = TRUE, direction = "<"), error = function(e) NULL)
  r2 <- tryCatch(pROC::roc(labels, s2, quiet = TRUE, direction = "<"), error = function(e) NULL)
  if (is.null(r1) || is.null(r2)) return(NA_real_)
  as.numeric(pROC::roc.test(r1, r2, method = "delong")$p.value)
}

# Define groups
obs <- test[[y_col]]  # ensure 'obs' defined for groups below
groups <- list(
  top_1  = list(labels = mk_labels_pct(obs, 0.01, "top"),    orient = "top"),
  bot_1  = list(labels = mk_labels_pct(obs, 0.01, "bottom"), orient = "bottom"),
  top_5  = list(labels = mk_labels_pct(obs, 0.05, "top"),    orient = "top"),
  bot_5  = list(labels = mk_labels_pct(obs, 0.05, "bottom"), orient = "bottom"),
  top_10 = list(labels = mk_labels_pct(obs, 0.10, "top"),    orient = "top"),
  bot_10 = list(labels = mk_labels_pct(obs, 0.10, "bottom"), orient = "bottom"),
  out_hi = list(labels = mk_labels_out(obs, "high"),         orient = "top"),
  out_lo = list(labels = mk_labels_out(obs, "low"),          orient = "bottom")
)
pretty_names <- c(
  top_1="Top 1%", bot_1="Bottom 1%",
  top_5="Top 5%", bot_5="Bottom 5%",
  top_10="Top 10%", bot_10="Bottom 10%",
  out_hi="Outlier High (>3)", out_lo="Outlier Low (<-3)"
)

# Compute AUCs & DeLong p-values (interaction vs cv / rv)
rows <- list(); comps <- list()
for (g in names(groups)) {
  lab <- groups[[g]]$labels
  if (length(unique(lab)) < 2L) next

  # Orientation: for bottom groups/outliers, invert scores so higher = positive
  s_cv  <- if (groups[[g]]$orient == "bottom") -pred_cv  else  pred_cv
  s_rv  <- if (groups[[g]]$orient == "bottom") -pred_rv  else  pred_rv
  s_int <- if (groups[[g]]$orient == "bottom") -pred_int else  pred_int

  a_cv  <- auc_ci(lab, s_cv)
  a_rv  <- auc_ci(lab, s_rv)
  a_int <- auc_ci(lab, s_int)

  p_int_vs_cv <- delong_p(lab, s_int, s_cv)
  p_int_vs_rv <- delong_p(lab, s_int, s_rv)

  rows[[length(rows)+1]] <- data.frame(group = pretty_names[g], model = "cvPRS",           AUC = a_cv[1],  CI_low = a_cv[2],  CI_high = a_cv[3])
  rows[[length(rows)+1]] <- data.frame(group = pretty_names[g], model = "rvPRS",           AUC = a_rv[1],  CI_low = a_rv[2],  CI_high = a_rv[3])
  rows[[length(rows)+1]] <- data.frame(group = pretty_names[g], model = "interaction.PRS", AUC = a_int[1], CI_low = a_int[2], CI_high = a_int[3])

  comps[[length(comps)+1]] <- data.frame(group = pretty_names[g], p_int_vs_cv = p_int_vs_cv, p_int_vs_rv = p_int_vs_rv)
}
res_auc  <- bind_rows(rows)
res_pval <- bind_rows(comps)

fwrite(res_auc,  file = "PRS_extremes_AUC_results.txt", sep = "\t")
fwrite(res_pval, file = "PRS_extremes_p_values.txt",   sep = "\t")

# Plot AUC (95% CI) with stars for interaction vs cvPRS
plot_df <- res_auc %>% left_join(res_pval, by = "group")

plot_df <- plot_df %>%
  mutate(sig_label = case_when(
    model == "interaction.PRS" & !is.na(p_int_vs_cv) & p_int_vs_cv < 5e-4 ~ "***",
    model == "interaction.PRS" & !is.na(p_int_vs_cv) & p_int_vs_cv < 5e-3 ~ "**",
    model == "interaction.PRS" & !is.na(p_int_vs_cv) & p_int_vs_cv < 5e-2 ~ "*",
    TRUE ~ ""
  ))

model_colors <- c("cvPRS" = "#1f77b4", "rvPRS" = "#d62728", "interaction.PRS" = "#9467bd")

order_groups <- c("Top 1%","Bottom 1%","Top 5%","Bottom 5%","Top 10%","Bottom 10%","Outlier High (>3)","Outlier Low (<-3)")
plot_df$group <- factor(plot_df$group, levels = order_groups)

stars_df <- plot_df %>%
  filter(model == "interaction.PRS", sig_label != "") %>%
  group_by(group) %>%
  summarise(sig_label = dplyr::first(sig_label),
            y_pos = max(CI_high, na.rm = TRUE) + 0.03,
            .groups = "drop")

p <- ggplot(plot_df, aes(x = group, y = AUC, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_color_manual(values = model_colors) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top",
        panel.grid.minor = element_blank()) +
  labs(x = "Group", y = "AUC (95% CI)", color = "Model") +
  geom_text(data = stars_df, inherit.aes = FALSE,
            aes(x = group, y = y_pos, label = sig_label),
            color = "black", size = 5, vjust = 0) +
  coord_cartesian(ylim = c(0.45, max(plot_df$CI_high, na.rm = TRUE) + 0.08))

ggsave("PRS_extremes_AUC_with_stars.pdf", p, width = 9.5, height = 5)

message("Done. Wrote: PRS_test_performance_metrics.txt / PRS_test_performance_table.pdf / ",
        "PRS_extremes_AUC_results.txt / PRS_extremes_p_values.txt / PRS_extremes_AUC_with_stars.pdf")
