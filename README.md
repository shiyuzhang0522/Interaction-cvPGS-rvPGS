# 1. Repository overview

This repository contains the analysis code accompanying the manuscript
**“Interactions between rare and common variant genetic risks in determining telomere length.”**
It provides the workflows and scripts required to reproduce the construction of common-variant and rare-variant polygenic scores (cvPGS and rvPGS), as well as downstream analyses.

---

# 2. Data preparation

**(1) Phenotype data (UK Biobank, approved applicants only)**

* Telomere length: Z-adjusted log-transformed T/S ratio (*UKB field ID: 22192*).
* Covariates: age, age², sex, and genetic principal components 1–10 (*UKB field ID: 22009*).

**(2) Genotype and sequencing data**

* Imputed genotype data (*Field ID: 22828; Imputation from genotype (WTCHG)*).
* Hard-call genotypes (*Field ID: 22418*).
* Whole-exome sequencing data (*Field ID: 23157; Population-level exome OQFE variants, pVCF format – final exome release*).
* These datasets were used to construct ensemble cvPGS and rvPGS. Quality control procedures are described in the *Supplementary Methods* of the manuscript.

**(3) Summary statistics**

* Common-variant GWAS summary statistics (generated with **REGENIE v4.1.gz**).
* Rare-variant association results (generated with **STAARpipeline v0.9.8**).
* Both sets of summary statistics are available via Figshare: [DOI: 10.6084/m9.figshare.29949551](https://doi.org/10.6084/m9.figshare.29949551).
* These files must be downloaded prior to running the analyses in this repository.

---

# 3. Code

**(1) RICE pipeline**

We followed the **RICE pipeline** to derive ensemble cvPGS and rvPGS, based on:

> Williams J, Chen T, Hua X, Wong W, Yu K, Kraft P, Li X, Zhang H.
> *Integrating common and rare variants improves polygenic risk prediction across diverse populations.*
> medRxiv 2024.11.05.24316779; [https://doi.org/10.1101/2024.11.05.24316779](https://doi.org/10.1101/2024.11.05.24316779)

For additional details, see the official repository: [RareVariantPRS](https://github.com/jwilliams10/RareVariantPRS.git).

**(2) Manuscript analyses**

Scripts used to generate figures and tables presented in the manuscript are included in this repository.

---

# 4. Usage

(1) RICE-pipeline/01_pheno_and_splits.R

This script processes phenotype data, removes missing telomere length samples, and generates **10-fold cross-validation splits** (train/tune/val), along with phenotype and covariate files for REGENIE.

### Example command

```bash
Rscript 01_pheno_and_splits.R
```

(2) RICE-pipeline/RICE-CV/  
This folder contains the scripts for constructing and evaluating the **ensemble cvPGS** using a cross-validation (CV) framework.  

Scripts are numbered (02–09) to reflect the recommended execution order.  
Each script takes the CV fold (and chromosome index where applicable) as input.  
Outputs from earlier steps serve as inputs for subsequent steps, culminating in the ensemble cvPGS for each CV fold.  

(3) RICE-pipeline/RICE-RV/  
For rare-variant association analyses, we recommend using the **STAARpipeline**.  
A step-by-step tutorial is available here: [STAARpipeline-Tutorial](https://github.com/xihaoli/STAARpipeline-Tutorial).  

**References**  
>[1] *A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies.* Nature Methods, 19(12), 1599–1611 (2022). PMID: 36303018; PMCID: PMC10008172; DOI: [10.1038/s41592-022-01640-x](https://doi.org/10.1038/s41592-022-01640-x).  
>[2] *Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale.* Nature Genetics, 52(9), >969–983 (2020). PMID: 32839606; PMCID: PMC7483769; DOI: [10.1038/s41588-020-0676-4](https://doi.org/10.>1038/s41588-020-0676-4).  

**Note**  
After running STAARpipeline and preparing gene-based rare variant masks, you can calculate ensemble rvPGS across cross-validation folds.

(4) Analysis in the Manuscript  
This folder contains scripts **01–10**, which reproduce the main figures and results presented in the manuscript.  

- **01** Compare the distributions, effect sizes, and model fits of cvPGS and rvPGS.  
- **02** Explore complementary features of cvPGS and rvPGS.  
- **03** Analyze correlations between scaled_cvPGS and scaled_rvPGS, including extreme groups.  
- **04** Create Miami plots combining GWAS and rare-variant collapsing results to show convergence of signals.  
- **05** Genome-wide interaction analysis between cvPGS and rvPGS, including group-based cutoffs.  
- **06** Test interactions between rare-variant burden masks and cvPGS, with visualization of slopes.  
- **07** Fine-mapping sentinel variants using [fine-mapping-inf](https://github.com/FinucaneLab/fine-mapping-inf) with region-specific LD and summary statistics.  
>  Ref: Cui R, Elzur RA, Kanai M, Ulirsch JC, Weissbrod O, Daly MJ, Neale BM, Fan Z, Finucane HK.  
>  *Improving fine-mapping by modeling infinitesimal effects.* Nat Genet. 2024;56(1):162–169. doi: [10.1038/s41588-023-01597-3](https://doi.org/10.1038/s41588-023-01597-3) 

- **08** Variant-level interaction analysis between common causal variants and rare-variant burden.  
- **09** Predict effects of **rs2853677** on transcription factor binding motifs using [motifbreakR](https://bioconductor.org/packages/motifbreakR)  
  (see [motifbreakR vignette](https://bioconductor.org/packages/devel/bioc/vignettes/motifbreakR/inst/doc/motifbreakR-vignette.html)).  
- **10** Combine cvPGS and rvPGS into additive and interaction models to assess predictive improvement.  

---

📬 **Contact**

For questions regarding this repository, please contact: **[shiyuzhang0522@gmail.com](mailto:shiyuzhang0522@gmail.com)**
