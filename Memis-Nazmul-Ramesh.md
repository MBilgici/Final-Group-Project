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
