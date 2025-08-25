#### 08. Variant-level Interaction between rare-variant burden per mask and causal variants after finemapping #############
################ Rscript: Common-variant × Rare-burden interactions (phenotype) ################
################ Author: Shelley ################################################################
################ Date: 2025-08-16 ################################################################

# Working directory
setwd("/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Interaction/Variant-level")

# ---- Packages -------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ---- Inputs ----------------------------------------------------------------------------------
PHENO_FILE   <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/merged_all_data.txt"
BURDEN_FILE  <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Interaction/Gene-level/ALL_validation_mask_burdens.txt"
GENO_FILE    <- "/Users/shelleyz/Downloads/Project-LTL-ScienceBulletin/Major.Revision.0724/ReAnalysis/Interaction/Variant-level/ALTdosage.IID_genotype_combined.tsv"
PHENO_COL    <- "scaled_resid_TS_log_Zadj"

# thresholds / housekeeping
MIN_TOTAL_N  <- 1000   # minimal complete-case N per variant×mask
MIN_B1_N     <- 50     # minimal number of burden carriers per mask
REPORT_EVERY <- 1000   # progress counter for scanning

# ---- Read data -------------------------------------------------------------------------------
dt <- fread(PHENO_FILE)  # expects: IID, scaled_resid_TS_log_Zadj, ...
burden_dt <- fread(BURDEN_FILE)        # IID + mask columns (counts)
genotype_dt <- fread(GENO_FILE)        # IID + variant genotype (0/1/2 ALT dosage)

# Keep genotypes for individuals in phenotype table
genotype_dt <- genotype_dt[IID %in% dt$IID]

# ---- Genotype: validate (0/1/2/NA) and build categorical calls --------------------------------
geno_cols <- setdiff(names(genotype_dt), "IID")

geno_calls012 <- copy(genotype_dt)
for (cc in geno_cols) {
  v <- suppressWarnings(as.integer(geno_calls012[[cc]]))
  bad <- !is.na(v) & !v %in% c(0L, 1L, 2L)
  if (any(bad)) {
    stop("Genotype column ", cc, " contains non-0/1/2/NA values. Examples: ",
         paste(unique(genotype_dt[[cc]][bad])[1:min(5, sum(bad))], collapse = ", "))
  }
  set(geno_calls012, j = cc, value = v)
}

map_levels <- function(v012) {
  factor(v012, levels = c(0L, 1L, 2L),
         labels = c("ref_hom", "het", "alt_hom"),
         ordered = TRUE)
}
geno_calls_cat <- copy(geno_calls012)
for (cc in geno_cols) {
  set(geno_calls_cat, j = cc, value = map_levels(geno_calls012[[cc]]))
}

# ---- Burden: convert to binary carrier (0/1), validate ---------------------------------------
burden_binary <- copy(burden_dt)
mask_cols <- setdiff(names(burden_binary), "IID")

# If counts are provided, convert to carriers (>=1 → 1; else 0; keep NA)
burden_binary[, (mask_cols) := lapply(.SD, function(x) {
  if (is.numeric(x)) as.integer(ifelse(is.na(x), NA, x > 0)) else as.integer(x)
}), .SDcols = mask_cols]

# Validate 0/1/NA
for (cc in mask_cols) {
  v <- suppressWarnings(as.integer(burden_binary[[cc]]))
  bad <- !is.na(v) & !(v %in% c(0L, 1L))
  if (any(bad)) {
    stop("Mask column ", cc, " contains non-binary values. Examples: ",
         paste(unique(burden_binary[[cc]][bad])[1:min(5, sum(bad))], collapse = ", "))
  }
  set(burden_binary, j = cc, value = v)
}

# Optional labeled version
burden_labeled <- copy(burden_binary)
lab <- function(x) factor(x, levels = c(0L, 1L), labels = c("non_carrier", "carrier"))
for (cc in mask_cols) set(burden_labeled, j = cc, value = lab(burden_binary[[cc]]))

# Per-mask summary
mask_summary <- rbindlist(lapply(mask_cols, function(cn) {
  v <- burden_binary[[cn]]
  n_missing <- sum(is.na(v))
  n_car     <- sum(v == 1L, na.rm = TRUE)
  n_non     <- sum(v == 0L, na.rm = TRUE)
  n_obs     <- n_car + n_non
  data.table(
    mask = cn, n_obs = n_obs, n_missing = n_missing,
    n_noncarrier = n_non, n_carrier = n_car,
    carrier_rate = ifelse(n_obs > 0, n_car / n_obs, NA_real_)
  )
}), use.names = TRUE)
setorder(mask_summary, -n_carrier)

# ---- Build analysis frame (Xall): IID + Y + genotypes (categorical) + burden (binary) ---------
stopifnot(PHENO_COL %in% names(dt))
Xall <- merge(dt[, .(IID, Y = get(PHENO_COL))], geno_calls_cat, by = "IID", all.x = FALSE, all.y = FALSE)
Xall <- merge(Xall, burden_binary, by = "IID", all.x = FALSE, all.y = FALSE)

# ---- Helper: fit a single variant×mask pair ---------------------------------------------------
fit_one_pair <- function(df, gcol, bcol,
                         std_levels = c("ref_hom","het","alt_hom"),
                         min_total_n = MIN_TOTAL_N,
                         min_b1_n = MIN_B1_N) {
  G <- df[[gcol]]
  B <- df[[bcol]]
  Y <- df[["Y"]]
  sub <- data.frame(Y = Y, G = G, B = B)
  sub <- sub[complete.cases(sub), , drop = FALSE]
  n <- nrow(sub)
  if (n < min_total_n) return(NULL)
  if (sum(sub$B == 1L) < min_b1_n) return(NULL)

  sub$G <- factor(sub$G, levels = std_levels)
  sub$G <- droplevels(sub$G)

  tb <- with(sub, table(G, B))  # 3×2 counts (may be smaller if levels absent)

  get_n <- function(g, b) if (g %in% rownames(tb) && b %in% colnames(tb)) tb[g, b] else 0L
  n_ref_B0 <- get_n("ref_hom", "0"); n_ref_B1 <- get_n("ref_hom", "1")
  n_het_B0 <- get_n("het",     "0"); n_het_B1 <- get_n("het",     "1")
  n_alt_B0 <- get_n("alt_hom", "0"); n_alt_B1 <- get_n("alt_hom", "1")

  n_ref <- n_ref_B0 + n_ref_B1
  n_het <- n_het_B0 + n_het_B1
  n_alt <- n_alt_B0 + n_alt_B1
  n_car <- n_ref_B1 + n_het_B1 + n_alt_B1

  fit_red  <- tryCatch(lm(Y ~ G + B, data = sub), error = function(e) NULL)
  if (is.null(fit_red)) return(NULL)
  fit_full <- tryCatch(lm(Y ~ G * B, data = sub), error = function(e) NULL)
  if (is.null(fit_full)) return(NULL)

  aov_cmp <- tryCatch(anova(fit_red, fit_full), error = function(e) NULL)
  p_int_2df <- if (!is.null(aov_cmp)) aov_cmp$`Pr(>F)`[2] else NA_real_

  co <- summary(fit_full)$coefficients
  rn <- rownames(co)
  get_term <- function(term) {
    if (term %in% rn) {
      c(beta = unname(co[term, "Estimate"]),
        se   = unname(co[term, "Std. Error"]),
        p    = unname(co[term, "Pr(>|t|)"]))
    } else {
      c(beta = NA_real_, se = NA_real_, p = NA_real_)
    }
  }

  t_Ghet   <- get_term("Ghet")
  t_Galt   <- get_term("Galt_hom")
  t_B      <- get_term("B")
  t_Ghet_B <- get_term("Ghet:B")
  t_Galt_B <- get_term("Galt_hom:B")

  data.table(
    variant     = gcol,
    mask        = bcol,
    n           = n,
    n_carrier   = n_car,
    n_ref_hom   = n_ref,
    n_het       = n_het,
    n_alt_hom   = n_alt,
    n_ref_B0    = as.integer(n_ref_B0),
    n_ref_B1    = as.integer(n_ref_B1),
    n_het_B0    = as.integer(n_het_B0),
    n_het_B1    = as.integer(n_het_B1),
    n_alt_B0    = as.integer(n_alt_B0),
    n_alt_B1    = as.integer(n_alt_B1),
    p_int_2df   = p_int_2df,
    beta_Ghet   = t_Ghet["beta"],   se_Ghet = t_Ghet["se"],   p_Ghet = t_Ghet["p"],
    beta_Galt   = t_Galt["beta"],   se_Galt = t_Galt["se"],   p_Galt = t_Galt["p"],
    beta_B      = t_B["beta"],      se_B    = t_B["se"],      p_B    = t_B["p"],
    beta_Ghet_B = t_Ghet_B["beta"], se_Ghet_B = t_Ghet_B["se"], p_Ghet_B = t_Ghet_B["p"],
    beta_Galt_B = t_Galt_B["beta"], se_Galt_B = t_Galt_B["se"], p_Galt_B = t_Galt_B["p"]
  )
}

# ---- Scan all variant×mask pairs --------------------------------------------------------------
total_tests <- length(geno_cols) * length(mask_cols)
cat("Planned tests:", total_tests, "\n")

res_list <- vector("list", total_tests)
k <- 0L; kept <- 0L
t_start <- proc.time()[3]

for (g in geno_cols) {
  Gv <- Xall[[g]]
  if (length(na.omit(unique(Gv))) < 2L) next  # need ≥2 genotype levels
  for (b in mask_cols) {
    k <- k + 1L
    out <- fit_one_pair(Xall, g, b)
    if (!is.null(out)) {
      kept <- kept + 1L
      res_list[[kept]] <- out
    }
    if (k %% REPORT_EVERY == 0L) {
      elapsed <- proc.time()[3] - t_start
      cat(sprintf("  tested %d / %d (kept %d) — elapsed %.1fs\n", k, total_tests, kept, elapsed))
      flush.console()
    }
  }
}
res_dt <- if (kept > 0L) rbindlist(res_list[seq_len(kept)], use.names = TRUE, fill = TRUE) else data.table()

# ---- Multiple testing (BH) --------------------------------------------------------------------
pcols <- c("p_int_2df", "p_Ghet", "p_Galt", "p_B", "p_Ghet_B", "p_Galt_B")
add_q <- function(DT, pcol) {
  if (pcol %in% names(DT)) {
    DT[, (paste0("q_", pcol)) := p.adjust(as.numeric(get(pcol)), method = "BH")]
  }
}
invisible(lapply(pcols, function(pc) add_q(res_dt, pc)))

# ---- Filter “most likely” interactions --------------------------------------------------------
need <- c("q_p_int_2df", "p_B", "p_Ghet", "p_Galt")
miss <- setdiff(need, names(res_dt))
if (length(miss)) stop("Missing required columns in res_dt: ", paste(miss, collapse = ", "))

hits <- res_dt[
  !is.na(q_p_int_2df) & q_p_int_2df < 0.05 &
  !is.na(p_B)    & p_B    < 2.5e-6  &
  !is.na(p_Ghet) & p_Ghet < 5e-8    &
  !is.na(p_Galt) & p_Galt < 5e-8
]
# Exclude HBB masks (measurement artifact concern)
hits <- hits[!grepl("HBB", mask, fixed = TRUE)]
setorder(hits, q_p_int_2df, p_int_2df, na.last = TRUE)

# ---- Save outputs -----------------------------------------------------------------------------
fwrite(res_dt, "scan_variantXburden_interactions.full.with_qvals.tsv", sep = "\t")
saveRDS(res_dt, "scan_variantXburden_interactions.full.with_qvals.rds")

fwrite(hits, "TERT_scan_variantXburden_interactions.filtered_hits.tsv", sep = "\t")
saveRDS(hits, "TERT_scan_variantXburden_interactions.filtered_hits.rds")

cat("Saved results:\n")
cat(" - Full results with q-values:  scan_variantXburden_interactions.full.with_qvals.tsv/.rds\n")
cat(" - Filtered hits:               TERT_scan_variantXburden_interactions.filtered_hits.tsv/.rds\n")
cat("   (Excluded masks containing 'HBB')\n")

if (nrow(hits)) {
  cat("Top filtered hits preview:\n")
  print(head(hits[, .(variant, mask,
                      n, n_carrier,
                      p_int_2df, q_p_int_2df,
                      beta_B, p_B,
                      beta_Ghet, p_Ghet,
                      beta_Galt, p_Galt,
                      beta_Ghet_B, p_Ghet_B,
                      beta_Galt_B, p_Galt_B)], 10))
} else {
  cat("No hits passed the strict filtering criteria.\n")
}

# ---- Visualization helpers --------------------------------------------------------------------
col_noncar <- "#1f77b4"  # blue
col_car    <- "#d62728"  # red

build_pair_df <- function(variant_col, mask_col) {
  df <- merge(dt[, .(IID, Y = get(PHENO_COL))],
              geno_calls_cat[, .(IID, G = get(variant_col))], by = "IID")
  df <- merge(df, burden_binary[, .(IID, B = get(mask_col))], by = "IID")
  df <- df[complete.cases(df[, .(Y, G, B)])]
  df[, G := factor(as.character(G), levels = c("ref_hom","het","alt_hom"))]
  df[, B := factor(as.integer(B), levels = c(0,1),
                   labels = c("non-carrier","carrier"))]
  df[]
}

summarise_cells <- function(d) {
  d[, .(
    n    = .N,
    mean = mean(Y),
    se   = sd(Y)/sqrt(.N),
    lwr  = mean(Y) - 1.96*sd(Y)/sqrt(.N),
    upr  = mean(Y) + 1.96*sd(Y)/sqrt(.N)
  ), by = .(G, B)]
}

pairwise_vs_ref <- function(d) {
  ref_vals <- d[G == "ref_hom" & B == "non-carrier", Y]
  d[, pval := t.test(Y, ref_vals, var.equal = FALSE)$p.value, by = .(G, B)]
  d[]
}

plot_interaction_tests <- function(variant_col, mask_col, results_table = NULL,
                                   out_file = NULL, title_prefix = "Interaction") {
  df   <- build_pair_df(variant_col, mask_col)
  if (!nrow(df)) return(NULL)

  cell <- summarise_cells(df)
  test_df <- pairwise_vs_ref(df)
  tests <- unique(test_df[, .(G, B, pval)])
  tests[, padj := p.adjust(pval, method = "BH")]
  cell <- merge(cell, tests, by = c("G", "B"))

  # global interaction q (if present)
  subt <- NULL
  if (!is.null(results_table) &&
      all(c("variant","mask","q_p_int_2df") %in% names(results_table))) {
    rr <- results_table[variant == variant_col & mask == mask_col]
    if (nrow(rr)) subt <- sprintf("FDR interaction q = %.2e", rr$q_p_int_2df[1])
  }

  p <- ggplot(cell, aes(x = G, y = mean, group = B, color = B)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.16, linewidth = 0.6) +
    scale_color_manual(values = c("non-carrier" = col_noncar, "carrier" = col_car),
                       name = "Burden") +
    labs(
      title = paste0(title_prefix, ": ", variant_col, " × ", mask_col),
      subtitle = subt,
      x = "Genotype", y = paste0("Mean ", PHENO_COL, " ± 95% CI")
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))

  # annotate adjusted p-values above each point
  p <- p + geom_text(aes(label = sprintf("q=%.2e", padj), y = upr + 0.05),
                     size = 3, vjust = 0, show.legend = FALSE)

  if (!is.null(out_file)) {
    ggsave(out_file, p, width = 6.5, height = 5.2)
  }
  p
}

# ---- Example: plot top three hits (if any) ----------------------------------------------------
p_list <- list()
if (nrow(hits) >= 1) p_list[[1]] <- plot_interaction_tests(hits$variant[1], hits$mask[1],
                                                           results_table = hits,
                                                           out_file = "interaction_top1_tests.png",
                                                           title_prefix = "Top hit")
if (nrow(hits) >= 2) p_list[[2]] <- plot_interaction_tests(hits$variant[2], hits$mask[2],
                                                           results_table = hits,
                                                           out_file = "interaction_top2_tests.png",
                                                           title_prefix = "Second hit")
if (nrow(hits) >= 3) p_list[[3]] <- plot_interaction_tests(hits$variant[3], hits$mask[3],
                                                           results_table = hits,
                                                           out_file = "interaction_top3_tests.png",
                                                           title_prefix = "Third hit")

p_list <- Filter(Negate(is.null), p_list)
if (length(p_list) > 0) {
  combined <- wrap_plots(p_list, nrow = 1, ncol = length(p_list), guides = "collect") &
    theme(legend.position = "bottom")
  ggsave("interaction_top3_meanCI_patchwork.pdf", combined, width = 16, height = 6)
}
