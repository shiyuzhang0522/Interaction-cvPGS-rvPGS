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

📬 **Contact**
For questions regarding this repository, please contact: **[shiyuzhang0522@gmail.com](mailto:shiyuzhang0522@gmail.com)**
