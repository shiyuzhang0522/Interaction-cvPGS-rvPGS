#### 06. Gene-level Interaction between rare-variant burden per mask and cvPGS #############
############ Author: Shelley ####################################################
############ Date: 2025-08-11 ###################################################

# --- Setup ---------------------------------------------------------------------
setwd("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Interaction/Gene-level")

library(data.table)
library(ggplot2)

# --- 1. Input data -------------------------------------------------------------
burden_dt <- fread("ALL_validation_mask_burdens.txt")   # IID + mask columns
dt        <- fread("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt")

# --- 2. Convert burden to binary carrier status --------------------------------
burden_binary <- copy(burden_dt)
mask_cols <- setdiff(names(burden_binary), "IID")
burden_binary[, (mask_cols) := lapply(.SD, function(x) as.integer(x > 0)), .SDcols = mask_cols]

mask_summary <- data.table(
  mask_id        = mask_cols,
  n_carriers     = sapply(burden_binary[, ..mask_cols], sum, na.rm = TRUE),
  n_non_carriers = nrow(burden_binary) - sapply(burden_binary[, ..mask_cols], sum, na.rm = TRUE)
)

fwrite(burden_binary, "ALL_validation_mask_burdens_binary.txt", sep="\t", quote=FALSE)
fwrite(mask_summary, "mask_carrier_summary.txt", sep="\t", quote=FALSE)

cat("Binary carrier matrix saved.\n")
cat("Mask carrier summary saved.\n")

# --- 3. Test interactions with cvPGS -------------------------------------------
burden_binary[, IID := as.character(IID)]
dt[, IID := as.character(IID)]
m <- merge(dt, burden_binary, by="IID", all=FALSE)

res_list <- lapply(mask_cols, function(mc) {
  x  <- m[[mc]]
  cv <- m[["scaled_cvPRS"]]
  y  <- m[["scaled_resid_TS_log_Zadj"]]
  
  keep <- complete.cases(x, cv, y)
  x <- x[keep]; cv <- cv[keep]; y <- y[keep]
  
  n_car <- sum(x == 1); n_non <- sum(x == 0)
  if (n_car == 0L || n_non == 0L) {
    return(data.table(mask_id=mc, n=length(x), n_carriers=n_car, n_noncarriers=n_non,
                      beta_cv=NA, se_cv=NA, t_cv=NA, p_cv=NA,
                      beta_mask=NA, se_mask=NA, t_mask=NA, p_mask=NA,
                      beta_int=NA, se_int=NA, t_int=NA, p_int=NA, r2_full=NA))
  }
  
  fit <- lm(y ~ cv * x)
  sm  <- summary(fit)$coefficients
  rn  <- rownames(sm)
  
  g <- function(name) if (name %in% rn) sm[name, ] else c(NA, NA, NA, NA)
  int_row <- if ("cv:x" %in% rn) "cv:x" else if ("x:cv" %in% rn) "x:cv" else NA
  
  b_cv   <- g("cv")
  b_mask <- g("x")
  b_int  <- if (!is.na(int_row)) sm[int_row, ] else c(NA, NA, NA, NA)
  
  data.table(mask_id=mc, n=length(x), n_carriers=n_car, n_noncarriers=n_non,
             beta_cv=b_cv[1], se_cv=b_cv[2], t_cv=b_cv[3], p_cv=b_cv[4],
             beta_mask=b_mask[1], se_mask=b_mask[2], t_mask=b_mask[3], p_mask=b_mask[4],
             beta_int=b_int[1], se_int=b_int[2], t_int=b_int[3], p_int=b_int[4],
             r2_full=summary(fit)$r.squared)
})

results <- rbindlist(res_list, fill=TRUE)
results[, p_int_BH := ifelse(is.na(p_int), NA_real_, p.adjust(p_int, method="BH"))]

fwrite(results, "rare.variant.burden_cvPGS_interactions_results.txt", sep="\t", quote=FALSE)
cat("Interaction results saved.\n")

# --- 4. Visualization of significant interactions ------------------------------
fit_slope <- function(y, x) {
  keep <- complete.cases(y, x); y <- y[keep]; x <- x[keep]
  if (length(y) < 3 || length(unique(x)) < 2)
    return(list(slope=NA, se=NA, n=length(y), lo=NA, hi=NA))
  sm <- summary(lm(y ~ x))$coefficients
  slope <- sm["x","Estimate"]; se <- sm["x","Std. Error"]; df <- length(y) - 2
  tcrit <- qt(0.975, df); lo <- slope - tcrit*se; hi <- slope + tcrit*se
  list(slope=slope, se=se, n=length(y), lo=lo, hi=hi)
}

sig_res <- results[!is.na(p_int_BH) & p_int_BH < 0.05]
stopifnot(nrow(sig_res) > 0)
ord_levels <- sig_res$mask_id[order(sig_res$p_int_BH)]

sig_masks <- as.character(unique(ord_levels))
m <- merge(dt[, .(IID, scaled_cvPRS, scaled_resid_TS_log_Zadj)],
           burden_binary[, c("IID", sig_masks), with=FALSE], by="IID")

rows <- lapply(sig_masks, function(msk) {
  xmask <- m[[msk]]
  y <- m$scaled_resid_TS_log_Zadj
  cv <- m$scaled_cvPRS
  nc <- fit_slope(y[xmask==0], cv[xmask==0])
  ca <- fit_slope(y[xmask==1], cv[xmask==1])
  data.table(mask_id=msk,
             group=c("non-carrier","carrier"),
             slope=c(nc$slope, ca$slope), se=c(nc$se, ca$se),
             n=c(nc$n, ca$n), ci_lo=c(nc$lo, ca$lo), ci_hi=c(nc$hi, ca$hi))
})
slope_long <- rbindlist(rows, use.names=TRUE)
slope_long[, mask_id := factor(mask_id, levels=ord_levels)]
slope_long[, group := factor(group, levels=c("non-carrier","carrier"))]

# p-value labels
lab_dt <- slope_long[, .(ymax=max(ci_hi,na.rm=TRUE),
                         yrng=diff(range(c(ci_lo,ci_hi),na.rm=TRUE))),
                     by=mask_id]
lab_dt <- merge(lab_dt, sig_res[, .(mask_id,p_int_BH)], by="mask_id", all.x=TRUE)
lab_dt[!is.finite(yrng)|is.na(yrng), yrng:=0]
lab_dt[, `:=`(y_lab=ymax+0.05*yrng,
              label=paste0("BH p = ", format.pval(p_int_BH, digits=2)))]
lab_dt[, mask_id := factor(mask_id, levels=ord_levels)]

p <- ggplot(slope_long, aes(x=group, y=slope)) +
  geom_line(aes(group=mask_id), color="grey60", linetype="dashed", linewidth=0.6) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi, color=group), width=0.15) +
  geom_point(aes(color=group), size=2.6) +
  geom_text(data=lab_dt, aes(x=1.5, y=y_lab, label=label),
            inherit.aes=FALSE, size=3.3) +
  facet_wrap(~ mask_id, ncol=4, scales="free_y") +
  scale_color_manual(values=c("non-carrier"="#2b8cbe", "carrier"="#de2d26")) +
  labs(x="Carrier status", y="Slope of cvPRS → outcome (95% CI)",
       title="cvPRS slopes by carrier status (within-group fits)",
       subtitle="Facets ordered by BH-adjusted p-value (most significant first)") +
  theme_bw(base_size=11) +
  theme(legend.position="top")

ggsave("predicted_slopes_within_group_CI_sigOrdered.pdf", p, width=12, height=10.5)
print(p)

fwrite(slope_long, "predicted.cvPGS.slopes_within_group_CI_sigOrdered_data.txt", sep="\t")
cat("Slope plot and data table saved.\n")
