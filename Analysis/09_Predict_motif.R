#### 09. prediction binding to motifs #############
################ Rscript: Predict TF motif impact (motifbreakR) ################
################ Author: Shelley ###############################################
################ Date: 2025-08-18 ##############################################

# Working dir
setwd("/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/Interaction/Sentinel.Variants/TF_change/breakmotifR")

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(motifbreakR)
  library(MotifDb)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh37)  # dbSNP build for GRCh37
  library(BSgenome.Hsapiens.UCSC.hg19)       # hg19 genome
  library(ggplot2)
  library(patchwork)
})

# ---- Inputs ------------------------------------------------------------------
target_rsid <- "rs2853677"
snps_rdata  <- sprintf("%s_snps.mb.Rdata", target_rsid)  # cached SNP GRanges

# Motif set bundled with motifbreakR (human PWMs)
data(motifbreakR_motif)  # object: motifbreakR_motif

# ---- Obtain SNP GRanges (load cache or build from rsID) ----------------------
if (file.exists(snps_rdata)) {
  load(snps_rdata)  # loads 'snps.mb'
  if (!exists("snps.mb")) stop("RData loaded but 'snps.mb' not found.")
} else {
  message("Building snps.mb from rsID via snps.from.rsid() ...")
  snps.mb <- snps.from.rsid(
    rsid          = target_rsid,
    dbSNP         = SNPlocs.Hsapiens.dbSNP155.GRCh37,
    search.genome = BSgenome.Hsapiens.UCSC.hg19
  )
  save(snps.mb, file = snps_rdata)
  message("Saved SNP GRanges to: ", snps_rdata)
}

# Quick check
stopifnot(length(snps.mb) >= 1)

# ---- Run motifbreakR ---------------------------------------------------------
# Use effect threshold at 0.85, log method (standard), scan all motifs
discover.motif <- motifbreakR(
  snpList  = snps.mb[1],           # single sentinel variant
  filterp  = FALSE,                # keep all p; threshold applies to effect size
  pwmList  = motifbreakR_motif,
  threshold= 0.85,
  method   = "log",
  verbose  = TRUE,
  BPPARAM  = BiocParallel::bpparam()
)

message("motifbreakR completed: ", length(discover.motif), " matches.")

# ---- Save results ------------------------------------------------------------
save(discover.motif, file = "breakmotifR.discover_motif_results.RData")

discover.motif.df <- as.data.frame(discover.motif, row.names = NULL)
fwrite(
  discover.motif.df,
  file = "breakmotifR.discover_motif_results.txt",
  sep = "\t", row.names = FALSE, quote = FALSE
)

# ---- Focused plots for selected TFs (effect == 'strong') ---------------------
genes_of_interest <- c("SMAD3", "RAD21", "SMC3", "CTCF")

# Helper: filter by TF geneSymbol and (optionally) effect
keep_genes <- function(gr, genes, effect = c("strong","weak")) {
  effect <- match.arg(effect)
  sym <- mcols(gr)$geneSymbol
  ok_gene <- !is.na(sym) & sym %in% genes
  if (effect == "any") {
    gr[ok_gene]
  } else {
    eff <- mcols(gr)$effect
    gr[ok_gene & !is.na(eff) & eff == effect]
  }
}

# Build per-gene GRanges (strong effects only)
per_gene <- lapply(genes_of_interest, function(g) keep_genes(discover.motif, g, effect = "strong"))
names(per_gene) <- genes_of_interest

# Plot and save each gene (PDF)
for (g in genes_of_interest) {
  gr <- per_gene[[g]]
  if (length(gr) > 0L) {
    pdf(sprintf("plot_%s_%s_strong.pdf", target_rsid, g), width = 6, height = 5)
    plotMB(results = gr, rsid = target_rsid, effect = "strong")
    dev.off()
  } else {
    message(sprintf("No 'strong' matches for %s.", g))
  }
}

# ---- Save per-gene subsets (RData + TXT if non-empty) ------------------------
for (g in genes_of_interest) {
  gr <- per_gene[[g]]
  obj_name <- sprintf("mb_%s_%s", target_rsid, g)
  assign(obj_name, gr)
  save(list = obj_name, file = sprintf("%s.Rdata", obj_name))
  if (length(gr) > 0L) {
    df <- as.data.frame(gr, row.names = NULL)
    fwrite(df, file = sprintf("%s.txt", obj_name), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

message("All done.")
