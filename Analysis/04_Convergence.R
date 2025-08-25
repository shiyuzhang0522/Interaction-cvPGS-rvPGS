#### 04. Miami plots for convergence of GWAS and rare-variant signals ###########
############# Author: Shelley ###################################################
############# Date: 2025-08-09 ##################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# --- 1) Load GWAS results ------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Error: Please provide a CV_number as the first argument.")
}
CV_number <- as.character(args[1])
cat(sprintf("Extracting GWAS results for chr1–chr22 for CV %s at %s.\n",
            CV_number, Sys.time()))

step2_dir <- "/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/RICE_pipeline/RICE-CV/CV-GWAS/Step2"

gwas_list <- lapply(1:22, function(chr) {
  f <- file.path(step2_dir,
                 sprintf("regenie_step2_CV%s_chr%d_TS_log_Zadj.regenie", CV_number, chr))
  if (file.exists(f)) fread(f) else NULL
})
gwas_cv <- rbindlist(gwas_list, use.names = TRUE, fill = TRUE)

# --- 2) Load rare-variant collapsing results -----------------------------------
rv_path <- sprintf(
  "/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/RV_Conditional_Analysis/CV%s/coding_sig.csv",
  CV_number
)
rv_coding_results <- fread(rv_path)[, .(
  Gene = `Gene name`,
  Chr,
  Category,
  P = `STAAR-O`
)]

# --- 3) Gene coordinates (GRCh37) ---------------------------------------------
gtf_file <- "/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/Miami.Plots/Homo_sapiens.GRCh37.87.gtf"

gtf <- fread(cmd = paste("grep -w 'gene' ", gtf_file), sep = "\t", header = FALSE)
setnames(gtf, c("seqname","source","feature","start","end",
                "score","strand","frame","attribute"))
gtf[, Gene := sub('.*gene_name "([^"]+)".*', '\\1', attribute)]
gene_coords <- unique(gtf[, .(Gene, Chr = seqname, start, end)])
gene_coords <- gene_coords[grepl("^[0-9]+$", Chr) & as.integer(Chr) %in% 1:22]
gene_coords <- gene_coords[Gene %in% unique(rv_coding_results$Gene)]
gene_coords[, BP := floor((start + end) / 2)]

rv_coding_results[, Chr := as.character(Chr)]
gene_coords[, Chr := as.character(Chr)]
rv_coding_results_annot <- merge(
  rv_coding_results,
  gene_coords[, .(Gene, Chr, BP)],
  by = c("Gene","Chr"),
  all.x = TRUE
)

# --- 4) Prepare GWAS and rare-variant datasets --------------------------------
chr_col <- if ("CHROM" %in% names(gwas_cv)) "CHROM" else "CHR"
bp_col  <- if ("GENPOS" %in% names(gwas_cv)) "GENPOS" else "BP"

gwas_df <- gwas_cv[LOG10P >= 3, .(
  Chr     = as.character(get(chr_col)),
  BP      = as.numeric(get(bp_col)),
  P       = 10^(-LOG10P),
  Label   = NA_character_,
  Analysis = "GWAS"
)]
gwas_df[, Chr := gsub("^chr","",Chr)][Chr %in% as.character(1:22)]
gwas_df <- gwas_df[is.finite(P) & P > 0 & !is.na(BP)]

rv_df <- rv_coding_results_annot[, .(
  Chr     = as.character(Chr),
  BP      = as.numeric(BP),
  P       = as.numeric(P),
  Label   = as.character(Gene),
  Analysis = "Rare"
)]
rv_df[, Chr := gsub("^chr","",Chr)][Chr %in% as.character(1:22)]
rv_df <- rv_df[is.finite(P) & P > 0 & !is.na(BP)]

# --- 5) Combine and compute signed log10P -------------------------------------
miami_df <- rbindlist(list(gwas_df, rv_df), use.names = TRUE, fill = TRUE)
miami_df[, logp := -log10(P)]
miami_df[Analysis == "GWAS", logp := -logp]

# --- 6) Genomic positions ------------------------------------------------------
miami_df[, Chr := gsub("^chr","",as.character(Chr))]
miami_df <- miami_df[Chr %in% as.character(1:22)]
miami_df[, Chr_num := as.integer(Chr)]
setorder(miami_df, Chr_num, BP)

chr_spans <- miami_df[, .(max_bp = max(BP,na.rm=TRUE)), by = Chr]
chr_spans[, Chr_num := as.integer(Chr)]
setorder(chr_spans, Chr_num)
chr_spans[, chr_start := shift(cumsum(max_bp), fill = 0)]
chr_spans[, center := chr_start + max_bp/2]

miami_df <- chr_spans[miami_df, on=.(Chr), nomatch=0L]
miami_df[, BPcum := BP + chr_start]
axis_df <- chr_spans[, .(Chr, center)]

# --- 7) Truncate extreme P-values ---------------------------------------------
p_trunc <- 1e-20
miami_df[, P := pmax(P,p_trunc)]
miami_df[, logp := -log10(P)]
miami_df[Analysis=="GWAS", logp := -logp]

gwas_pts <- miami_df[Analysis=="GWAS"]
rare_pts <- miami_df[Analysis=="Rare"]

rare_thresh <- -log10(2.5e-6)
gwas_thresh <- -log10(5e-8)
ylim_max <- max(abs(range(miami_df$logp, na.rm=TRUE)),
                rare_thresh, gwas_thresh)

chr_colors <- ifelse(axis_df$Chr %in% as.character(seq(1,22,2)),
                     "#1F77B4","#D62728")

# --- 8) Plot base Miami -------------------------------------------------------
p <- ggplot() +
  geom_point(data=gwas_pts,
             aes(x=BPcum,y=logp,color=as.factor(Chr)),
             size=0.5,alpha=0.5,show.legend=FALSE) +
  geom_point(data=rare_pts,
             aes(x=BPcum,y=logp,color=as.factor(Chr)),
             size=1.0,alpha=0.85,show.legend=FALSE) +
  geom_hline(yintercept= rare_thresh, linetype="dashed", color="grey40") +
  geom_hline(yintercept=-gwas_thresh, linetype="dashed", color="grey40") +
  scale_x_continuous(breaks=axis_df$center, labels=axis_df$Chr) +
  scale_color_manual(values=chr_colors) +
  coord_cartesian(ylim=c(-ylim_max,ylim_max)) +
  labs(x="Chromosome", y=expression(-log[10](P))) +
  theme_bw(base_size=13) +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank())

# --- 9) Add convergence lines -------------------------------------------------
prune_window <- function(dt, kb=1e4) {
  if (!nrow(dt)) return(dt)
  out <- vector("list",0L); taken <- rep(FALSE,nrow(dt))
  for (i in seq_len(nrow(dt))) {
    if (taken[i]) next
    out[[length(out)+1L]] <- dt[i]
    taken <- taken | (dt$Chr==dt$Chr[i] & abs(dt$BP-dt$BP[i]) <= kb)
  }
  rbindlist(out)
}

rv_sig <- rare_pts[P <= 2.5e-6, .(Chr,BP,BPcum,P,logp,Label)]
setorder(rv_sig,Chr,P)
rv_leads <- prune_window(rv_sig, kb=1e4)

gwas_pool <- gwas_pts[P <= 5e-8, .(Chr,BP,BPcum,P,logp)]
if (!nrow(gwas_pool)) gwas_pool <- gwas_pts[P < 1e-3, .(Chr,BP,BPcum,P,logp)]
setorder(gwas_pool,Chr,P)
gwas_leads <- prune_window(gwas_pool, kb=1e4)

setkey(gwas_leads, Chr, BP)
nearest_gwas <- gwas_leads[rv_leads, roll="nearest", on=.(Chr,BP)]
pairs <- cbind(rv_leads[, .(Chr,BP_rv=BP,BPcum_rv=BPcum,Gene=Label)],
               nearest_gwas[, .(BP_gwas=BP,BPcum_gwas=BPcum,P_gwas=P)])
pairs[, dist_bp := abs(BP_rv-BP_gwas)]
pairs <- pairs[dist_bp <= 1e4]

lines_df <- pairs[, .(Chr, x=(BPcum_rv+BPcum_gwas)/2)]

p <- p + geom_vline(data=lines_df, aes(xintercept=x),
                    color="#6A5ACD", linetype="dotdash",
                    linewidth=0.7, alpha=0.6)

# --- 10) Save output ----------------------------------------------------------
out_dir  <- "/gpfs/hpc/home/lijc/zhangsy/Sci.Bull.Revision/Re-analysis/Miami.Plots"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
out_file <- file.path(out_dir,
                      sprintf("Miami_plot_CV%s_with_convergence.pdf", CV_number))
ggsave(out_file, p, width=11, height=6, device=cairo_pdf)
cat(sprintf("Plot with convergence lines saved to %s\n", out_file))
