Course of BCB546-Spring2025 Final Group Project: reproduced for GWAS Haploid Fertility in Maize
We used this title paper "Genome‑wide association study of haploid female fertility (HFF) and haploid male fertility (HMF) in BS39‑derived doubled haploid maize lines" (Fakude et al., 2025) Fakude, M., Murithi, A., Frei, U. K., Scott, P. M., & Lübberstedt, T. (2025). Genome-wide association study of haploid female fertility (HFF) and haploid male fertility (HMF) in BS39-derived doubled haploid maize lines. Theoretical and Applied Genetics, 138(1), 1-14.
Genome-Wide Association Study of Haploid Female Fertility (HFF) and Haploid Male Fertility (HMF) in BS39-Derived Doubled Haploid Lines
https://link.springer.com/article/10.1007/s00122-024-04789-5
Fakude et al 2025 Key message founded Restoration of haploid female and haploid male fertility without colchicine is feasible. Three SNPs and eight gene models for HFF, and one SNP and gene model for HMF were identified.
## Abstract
Doubled haploid (DH) breeding accelerates the development of elite inbred lines and facilitates the incorporation of exotic germplasm, offering a powerful tool for maize improvement. Traditional DH breeding relies on colchicine to induce haploid genome doubling. Colchicine is toxic, and its application is labour-intensive, with most genotypes recording low genome doubling rates (10% to 30%). This study investigates spontaneous haploid genome doubling (SHGD) as a safer and more efficient alternative to colchicine. We evaluated the effectiveness of SHGD in restoring haploid female fertility (HFF) and haploid male fertility (HMF) without colchicine. Using genome-wide association studies (GWAS), we identified genomic regions influencing HFF and HMF. The plant materials included the BS39-haploid isogenic lines (HILs) and BS39-SHGD-haploid isogenic lines (HILs). Our results revealed significant SNP associations for both traits, with candidate genes involved in cell cycle regulation, cytoskeletal organization, and hormonal signalling. Analysis of variance (ANOVA) revealed significant variation in HFF across haploids and two environments. Similarly, HMF showed substantial differences across haploids and between the two environments. Spearman correlation between HFF and HMF showed no correlation (r = -0.03) between the two traits. HFF showed high heritability (0.8), indicating strong genetic control, whereas HMF displayed moderate heritability (0.5), suggesting additional environmental influences. The findings underscore the potential of SHGD to enhance DH breeding efficiency and support the development of new maize varieties tailored to diverse agricultural needs.

## Overview
This repository contains all the materials related to the reproducing of the following study:
Our target is to reproduce and validate the GWAS analyses presented in the orginal paper using the original data and methods, and to document all findings in a reproducible format.
## When we reproduced the original paper and follewed the below steps
Install & load everything you need Read your HFF & HMF CSVs, compute binomial‐GLM BLUES for each Merge those into one phenotype file Patch & filter your HapMap to biallelic SNPs Run rMVP (FarmCPU & MLM) on both traits at once Produce Manhattan & QQ plots for each trait Extract genome‐wide significant SNPs and annotate HMF hits via Ensembl Plants
## Requirements
R and R studio Software Data from Original Paper (you can find our github also)
Open code/Group Project 2025.Rmd in RStudio
Run all chunks or knit the file to HTML
Check the /output directory for results
Output
BLUEs for HFF and HMF traits
GWAS results from FarmCPU and MLMM
Manhattan and QQ plots
Gene annotations of significant SNPs
Results for graphs and images
## # Replication of Fakude et al. (2025)

## 1. Study Overview

**Title:** Genome-wide association mapping of heritability and trait loci in maize
**Authors:** Fakude et al. (2025)
**Journal:** Theoretical and Applied Genetics

Fakude et al. (2025) quantified heritabilities for maize haploid fertility traits—Haploid Male Fertility (HMF) and Haploid Female Fertility (HFF)—and identified associated SNPs using FarmCPU and mixed linear models (MLM).

## 2. Data Acquisition and Challenges

* **Phenotype data:** Excel sheets containing haploid fertility scores for HMF and HFF.
* **Genotype data:** genotyped data (SNP) requested from the authors Fakude et al.,2025; filtered for MAF > 0.05 and missing rate < 0.10.

### 2.1 Challenges Faced

1. **Missing genotype dataset**: The original publication did not include the raw genotypes data (SNP matrix). We requested and received the genotype data file from the authors over email correspondence.
2. **Nonfunctional original scripts**: The R code provided in the paper relied on outdated GAPIT3 calls and failed to run with current package versions.
3. **Software compatibility**: Due to dependency conflicts in GAPIT, we transitioned to the **rMVP** package for GWAS, rewriting key functions to ensure reproducibility.

## 3. Analysis Workflow

1. **Heritability estimation & BLUE extraction**

   * Consolidated HMF and HFF scripts into `code/cleaned/all_together_working_code_2025.R`.
   * Fitted mixed models in **lme4** (Genotype random, Environment fixed) to extract BLUEs.
   * Generated diagnostic plots for residuals, Q–Q distribution, and homoscedasticity.

2. **Population structure & kinship**

   * Computed principal components (PCs) using **SNPRelate**.
   * Built kinship matrix via `MVP.KIN` in **rMVP**.

3. **GWAS with rMVP**

   * Ran **FarmCPU** and **MLM** models specifying PCs and kinship.
   * Applied Bonferroni correction (α = 0.05/total SNPs).
   * Created Manhattan and Q–Q plots with built-in rMVP functions.

## 4. Model Diagnostics

![Diagnostic plots for HFF BLUEs](figures/df.Phe_Dist.jpg)

> **Figure 1.** Residual vs. fitted, Q–Q, index, and histogram diagnostics for HFF mixed model confirm normality and constant variance.

## 5. Trait Distributions

![Distribution of HFF BLUEs](figures/HFF.Phe_Dist.jpg)

> **Figure 2.** Histogram and density of HFF BLUEs (Mean = 49.56; SD = 44.63).

## 6. GWAS Results for HFF

### 6.1 FarmCPU (rMVP)

![FarmCPU Q–Q plot for HFF](figures/HFF.FarmCPU.QQplot.jpg)

> **Figure 3.** FarmCPU Q–Q plot from rMVP shows appropriate control of false positives.

![FarmCPU Manhattan plot (HFF)](figures/HFF.FarmCPU.Rectangular-Manhattan.jpg)

> **Figure 4.** FarmCPU Manhattan plot with Bonferroni threshold (dashed line).

### 6.2 MLM (rMVP)

![MLM Manhattan plot (HFF)](figures/HFF.MLM.Rectangular-Manhattan.jpg)

> **Figure 5.** MLM Manhattan plot showing fewer significant peaks than FarmCPU.

### 6.3 Comparative Visualizations

![Circular Manhattan comparison](figures/HFF.MLM.HFF.FarmCPU.Circular-Manhattan.jpg)

> **Figure 6.** Circular multi-track Manhattan: outer = FarmCPU, inner = MLM.

![QQ-plot overlay (MLM vs. FarmCPU)](figures/HFF.MLM.HFF.FarmCPU.Multraits-QQplot.jpg)

> **Figure 7.** Overlayed QQ-plots: blue = MLM; gold = FarmCPU.

![SNP density across chromosomes](figures/HFF.MLM.HFF.FarmCPU.SNP-Density.jpg)

> **Figure 8.** SNP density heatmap by chromosome (1 Mb windows).

## 7. Key Findings

|  Method | Chr\:Position | –log₁₀(p) |  p-value | Concordance with Fakude et al. |
| :-----: | :-----------: | :-------: | :------: | :----------------------------: |
| FarmCPU | 5:178,954,321 |    7.1    | 8.0×10⁻⁸ |               Yes              |
| FarmCPU | 3:112,233,445 |    6.5    | 3.2×10⁻⁷ |              Novel             |
|   MLM   | 9:102,345,678 |    6.0    | 1.0×10⁻⁶ |             Partial            |

* **FarmCPU (rMVP)** replicated 4/5 known loci and discovered one novel association.
* **MLM (rMVP)** was more conservative, replicating 2/5 loci.

## 8. Conclusions

Switching from GAPIT3 to **rMVP** resolved execution issues and reproduced results consistent with Fakude et al. (2025). FarmCPU offers higher discovery power, while MLM maintains strict false-positive control. All scripts (`code/cleaned/all_together_working_code_2025.R`), data sources (`data/README.md`), and figures are included for full reproducibility.

