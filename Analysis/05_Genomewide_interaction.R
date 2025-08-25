#### 05. Genome-wide interaction analysis between cvPRS and rvPRS ################
############ Author: Shelley ####################################################
############ Date: 2025-08-11 ###################################################

# --- Setup ---------------------------------------------------------------------
library(data.table)
library(dplyr)
library(broom)
library(sandwich)
library(lmtest)
library(ggplot2)
library(purrr)
library(patchwork)

workdir <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Interaction/genome-wide"
infile  <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt"
setwd(workdir)

dt <- fread(infile)

# --- 1. Additive vs interaction model ------------------------------------------
logfile <- file.path(workdir, "interaction_vs_additive_Ftest.log")
capture.output({
  cat("### Comparing additive vs interaction model ###\n\n")
  
  fit_add <- lm(scaled_resid_TS_log_Zadj ~ scaled_cvPRS + scaled_rvPRS, data = dt)
  fit_int <- lm(scaled_resid_TS_log_Zadj ~ scaled_cvPRS * scaled_rvPRS, data = dt)
  
  a <- anova(fit_add, fit_int)
  print(a)
  
  r2_add   <- summary(fit_add)$r.squared
  r2_int   <- summary(fit_int)$r.squared
  delta_r2 <- r2_int - r2_add
  
  F_stat <- a$F[2]; p_val <- a$`Pr(>F)`[2]
  ss_diff <- a$`Sum of Sq`[2]
  rss_int <- deviance(fit_int)
  partial_r2 <- ss_diff / (ss_diff + rss_int)
  
  cat("\n### Model fit statistics ###\n")
  cat(sprintf("Additive R²:     %.6f\n", r2_add))
  cat(sprintf("Interaction R²:  %.6f\n", r2_int))
  cat(sprintf("ΔR² (int - add): %.6f\n", delta_r2))
  cat(sprintf("F = %.3f, p = %.3g\n", F_stat, p_val))
  cat(sprintf("Partial R² (interaction): %.6f\n", partial_r2))
}, file = logfile)

cat("Log saved to:", logfile, "\n")

summary(fit_int)

# --- 2. Visualize continuous interaction ---------------------------------------
b_cv  <- coef(fit_int)["scaled_cvPRS"]
b_rv  <- coef(fit_int)["scaled_rvPRS"]
b_int <- coef(fit_int)["scaled_cvPRS:scaled_rvPRS"]

rv_vals <- seq(-3, 3, 0.1)
cv_vals <- seq(-3, 3, 0.1)
slope_cv <- b_cv + b_int * rv_vals
slope_rv <- b_rv + b_int * cv_vals

p_slope_cv <- ggplot(data.frame(rvPRS = rv_vals, slope = slope_cv),
                     aes(rvPRS, slope)) +
  geom_line(size = 1.2, color = "#99000d") +
  geom_hline(yintercept = b_cv, linetype = "dashed") +
  labs(x = "rvPGS (SD units)", y = "Slope of cvPGS",
       title = "Change in cvPGS effect across rvPGS levels") +
  theme_minimal(base_size = 14)

p_slope_rv <- ggplot(data.frame(cvPRS = cv_vals, slope = slope_rv),
                     aes(cvPRS, slope)) +
  geom_line(size = 1.2, color = "#0c59cdff") +
  geom_hline(yintercept = b_rv, linetype = "dashed") +
  labs(x = "cvPGS (SD units)", y = "Slope of rvPGS",
       title = "Change in rvPGS effect across cvPGS levels") +
  theme_minimal(base_size = 14)

p_beta <- coef(summary(fit_int))["scaled_cvPRS:scaled_rvPRS", "Pr(>|t|)"]

combined_slope <- p_slope_cv + p_slope_rv +
  plot_layout(nrow = 2) +
  plot_annotation(
    subtitle = sprintf(
      "Interaction beta = %.4f; p = %s",
      b_int, format.pval(p_beta, digits = 3, eps = 1e-16)
    )
  )

ggsave("interaction_slopes_with_p_beta_vertical.pdf", combined_slope,
       width = 6, height = 12)

# --- 3. Binary×binary group-based interactions ---------------------------------
cutoffs <- c("highest1%", "highest5%", "highest10%", "+1SD", "+2SD", "+3SD",
             "lowest1%",  "lowest5%",  "lowest10%",  "-1SD", "-2SD", "-3SD")

get_threshold <- function(x, cut) {
  qs <- quantile(x, probs = c(.01,.05,.10,.90,.95,.99), na.rm = TRUE)
  switch(cut,
    "highest1%"  = list(thr=qs[["99%"]], dir="greater_or_equal"),
    "highest5%"  = list(thr=qs[["95%"]], dir="greater_or_equal"),
    "highest10%" = list(thr=qs[["90%"]], dir="greater_or_equal"),
    "+1SD"=list(thr=1, dir="greater_or_equal"),
    "+2SD"=list(thr=2, dir="greater_or_equal"),
    "+3SD"=list(thr=3, dir="greater_or_equal"),
    "lowest1%"   = list(thr=qs[["1%"]], dir="less_or_equal"),
    "lowest5%"   = list(thr=qs[["5%"]], dir="less_or_equal"),
    "lowest10%"  = list(thr=qs[["10%"]],dir="less_or_equal"),
    "-1SD"=list(thr=-1, dir="less_or_equal"),
    "-2SD"=list(thr=-2, dir="less_or_equal"),
    "-3SD"=list(thr=-3, dir="less_or_equal"),
    stop("Unknown cutoff: ", cut))
}
flag_group <- function(x, thr, dir) if (dir=="greater_or_equal") as.integer(x>=thr) else as.integer(x<=thr)

run_one_cutoff <- function(cut, dat) {
  th_cv <- get_threshold(dat$scaled_cvPRS, cut)
  th_rv <- get_threshold(dat$scaled_rvPRS, cut)
  df <- dat %>%
    mutate(cv_flag=flag_group(scaled_cvPRS,th_cv$thr,th_cv$dir),
           rv_flag=flag_group(scaled_rvPRS,th_rv$thr,th_rv$dir),
           cv_bin=factor(cv_flag,levels=c(0,1),labels=c("Non-group","Group")),
           rv_bin=factor(rv_flag,levels=c(0,1),labels=c("Non-group","Group")))
  fit <- lm(scaled_resid_TS_log_Zadj ~ cv_bin * rv_bin, data=df)
  tt <- broom::tidy(fit)
  b_cv <- tt %>% filter(term=="cv_binGroup") %>% slice_head(n=1)
  b_rv <- tt %>% filter(term=="rv_binGroup") %>% slice_head(n=1)
  b_int<- tt %>% filter(term=="cv_binGroup:rv_binGroup") %>% slice_head(n=1)
  tibble(cutoff=cut, beta_cv=b_cv$estimate, se_cv=b_cv$std.error,
         beta_rv=b_rv$estimate, se_rv=b_rv$std.error,
         beta_int=b_int$estimate,se_int=b_int$std.error,p_int=b_int$p.value)
}

res_tbl <- map_dfr(cutoffs, ~ run_one_cutoff(.x, dt)) %>%
  mutate(p_int_BH = p.adjust(p_int, method="BH"))
fwrite(res_tbl, file.path(workdir,"binary_interactions_specified_cutoffs.tsv"),
       sep="\t", quote=FALSE)

# --- 4. Forest plot of interaction betas ---------------------------------------
plot_order <- c("-2SD","-1SD","lowest10%","lowest5%","lowest1%",
                "+3SD","+2SD","+1SD","highest10%","highest5%","highest1%")

res_plot <- res_tbl %>%
  filter(cutoff %in% plot_order) %>%
  mutate(cutoff=factor(cutoff, levels=plot_order, ordered=TRUE),
         ci_lower=beta_int-1.96*se_int, ci_upper=beta_int+1.96*se_int,
         col_class=case_when(is.na(p_int_BH)|p_int_BH>0.5 ~ "gray",
                             beta_int>0 ~ "pos", beta_int<0 ~ "neg"),
         sig_lab=case_when(p_int_BH<1e-3~"***",p_int_BH<0.05~"**",
                           p_int_BH<0.10~"*",TRUE~""))

col_map <- c(pos="#de2d26", neg="#2b8cbe", gray="grey60")
span <- diff(range(c(res_plot$ci_lower,res_plot$ci_upper),na.rm=TRUE))

p <- ggplot(res_plot, aes(x=beta_int, y=cutoff, color=col_class)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey70") +
  geom_errorbarh(aes(xmin=ci_lower,xmax=ci_upper), height=0.22, linewidth=0.6) +
  geom_point(size=3) +
  geom_text(aes(label=sig_lab), nudge_x=0.03*span, vjust=0.35, size=4, na.rm=TRUE) +
  scale_color_manual(values=col_map, breaks=c("pos","neg","gray"),
                     labels=c("positive β_int","negative β_int","p > 0.5 or NA")) +
  labs(x=expression(beta[interaction]~"(95% CI)"), y="PGS cutoff", color=NULL,
       title="Interaction β across PGS cutoffs",
       subtitle="Red=positive; Blue=negative; Gray=BH p > 0.5. Stars: * <0.10, ** <0.05, *** <1e-3") +
  theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(), legend.position="right")

ggsave("forest.interaction_beta_by_sign_ordered_p05gray_no_disp.pdf", p,
       width=8, height=12)
