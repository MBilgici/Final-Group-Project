# Course of BCB546-Spring2025 Final Group Project: Replication for GWAS Haploid Fertility in Maize
We used this title paper "Genome‑wide association study of haploid female fertility (HFF) and haploid male fertility (HMF) in BS39‑derived doubled haploid maize lines" (Fakude et al., 2025)
Fakude, M., Murithi, A., Frei, U. K., Scott, P. M., & Lübberstedt, T. (2025). Genome-wide association study of haploid female fertility (HFF) and haploid male fertility (HMF) in BS39-derived doubled haploid maize lines. Theoretical and Applied Genetics, 138(1), 1-14.
## Overview
This repository contains all the materials related to the replication of the following study:

Our target is to replicate and validate the GWAS analyses presented in the orginal paper using the original data and methods, and to document all findings in a reproducible format.

## Repository Structure
1-README.md
2-Fakude et al., 2025.md (Summary and analysis of the original (Fakude et al., 2025) paper)
3-R code (R scripts and RMarkdown for data processing and analysis)
gwas_hff_hmf_analysis.Rmd
4-data
data sources
5-Our results output (Generated plots, summary tables, and annotations)
GWAS
Annotated_HFF_SNPs.csv
Annotated_HMF_SNPs.csv

## When we reproduced the original paper and follewed the below steps

### Requirements
R and R studio Software
Data from Original Paper (you can find our github also)

### Instructions
1. Clone this repository:
   ```bash
  git clone https://github.com/MBilgici/Final-Group-Project.git
  cd maize-haploid-gwas-replication
   ```
2. Open `code/gwas_hff_hmf_analysis.Rmd` in RStudio
3. Run all chunks or knit the file to HTML
4. Check the `/output` directory for results

### Output
- BLUEs for HFF and HMF traits
- GWAS results from FarmCPU and MLMM
- Manhattan and QQ plots
- Gene annotations of significant SNPs

