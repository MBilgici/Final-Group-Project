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
   ```{r libraries}
library(lme4)
library(ggResidpanel)
library(dplyr)
library(ggplot2)
library(emmeans)
library(arm)
library(tidyverse)
library(bestNormalize)
library(car)
library(qqman)
library(CMplot)
library(biomaRt)
```

## HMF Data Processing and BLUEs

```{r hmf_data}
hmf_raw <- read.csv("HMF_LAST.csv")
hmf_raw$PI <- as.factor(hmf_raw$PI)
hmf_raw$Rep <- as.factor(hmf_raw$Rep)

hmf_agg <- hmf_raw %>%
  group_by(PI, Rep) %>%
  summarize(Fertile = sum(Fertile), Sterile = sum(Sterile), .groups = "drop") %>%
  mutate(HMF = Fertile / (Fertile + Sterile))

min_nonzero <- min(hmf_agg$HMF[hmf_agg$HMF > 0], na.rm = TRUE) / 2
hmf_agg <- hmf_agg %>%
  mutate(HMF.mod = ifelse(HMF == 0, min_nonzero, HMF))

glm_model <- glm(cbind(Fertile, Sterile) ~ Rep + PI, data = hmf_agg, family = binomial)
glm_blues <- data.frame(emmeans(glm_model, "PI"))
glm_blues$HMF_BLUE <- arm::invlogit(glm_blues$emmean) * 100
```

## HFF Data Processing and BLUEs

```{r hff_data}
hff_data <- read.csv("HFF4.csv")
hff_data$PI <- as.factor(hff_data$PI)
hff_data$Rep <- as.factor(hff_data$Rep)
hff_data$HFF <- as.numeric(gsub("[^0-9.]", "", hff_data$HFF))
hff_data <- hff_data[!is.na(hff_data$HFF), ]

norm_result <- bestNormalize(hff_data$HFF)
hff_data$HFF_trans <- norm_result$chosen_transform$x.t

hff_model <- lm(HFF_trans ~ Rep + PI, data = hff_data)
hff_blues <- data.frame(emmeans(hff_model, "PI"))
hff_blues$HFF_BLUE <- hff_blues$emmean
```

## Summary Statistics

```{r summary_stats}
summary(hff_data$HFF)
sd(hff_data$HFF, na.rm = TRUE)
summary(hmf_agg$HMF)
sd(hmf_agg$HMF, na.rm = TRUE)
```

## GWAS Analysis with GAPIT

```{r gapit_gwas, eval=FALSE}
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

geno_data <- read.delim("genotype_data.txt", sep = "\t", header = TRUE)

geno_data <- geno_data[nchar(geno_data$alleles) == 3 & grepl("^[ACGT]/[ACGT]$", geno_data$alleles), ]
geno_data[, 12:ncol(geno_data)][geno_data[, 12:ncol(geno_data)] == "N"] <- NA

GAPIT(
  Y = glm_blues[, c("PI", "HMF_BLUE")],
  G = geno_data,
  PCA.total = 3,
  model = c("FarmCPU", "MLMM"),
  Multiple_analysis = TRUE,
  SNP.MAF = 0.01,
  file.output = TRUE
)

GAPIT(
  Y = hff_blues[, c("PI", "HFF_BLUE")],
  G = geno_data,
  PCA.total = 3,
  model = c("FarmCPU", "MLMM"),
  Multiple_analysis = TRUE,
  SNP.MAF = 0.01,
  file.output = TRUE
)
```

## Manhattan and QQ Plots

```{r manhattan_plots, fig.width=10, fig.height=6}
farmcpu_hff <- read.csv("GAPIT.Association.GWAS_Results.FarmCPU.HFF_BLUE.csv")
farmcpu_hmf <- read.csv("GAPIT.Association.GWAS_Results.FarmCPU.HMF_BLUE.csv")

par(mfrow = c(1, 2))
manhattan(farmcpu_hff, chr="Chr", bp="Pos", p="P.value", snp="SNP", main="FarmCPU - HFF")
manhattan(farmcpu_hmf, chr="Chr", bp="Pos", p="P.value", snp="SNP", main="FarmCPU - HMF")
```

```{r qq_plots, fig.width=10, fig.height=5}
par(mfrow = c(1, 2))
qq(farmcpu_hff$P.value, main = "QQ Plot - HFF")
qq(farmcpu_hmf$P.value, main = "QQ Plot - HMF")
```

## Top SNPs

```{r top_snps}
alpha <- 0.05 / nrow(farmcpu_hff)
top_snps_hff <- farmcpu_hff %>% filter(P.value <= alpha)
top_snps_hmf <- farmcpu_hmf %>% filter(P.value <= alpha)

head(top_snps_hff)
head(top_snps_hmf)
```

## Gene Annotations Using Ensembl Plants

```{r gene_annotation}
ensembl <- useMart("plants_mart", host = "https://plants.ensembl.org", dataset = "zmays_eg_gene")

annotate_snps <- function(snp_df, window = 200000) {
  gene_annotations <- list()

  for (i in 1:nrow(snp_df)) {
    chr <- as.character(snp_df$Chr[i])
    pos <- snp_df$Pos[i]

    genes <- getBM(
      attributes = c("chromosome_name", "start_position", "end_position", 
                     "ensembl_gene_id", "external_gene_name", "description", "gene_biotype"),
      filters = c("chromosome_name", "start", "end"),
      values = list(chr, pos - window, pos + window),
      mart = ensembl
    )

    if (nrow(genes) > 0) {
      genes$SNP <- snp_df$SNP[i]
      genes$SNP_Position <- pos
      gene_annotations[[i]] <- genes
    }
  }

  do.call(rbind, gene_annotations)
}

# Run on top SNPs
annotated_hff <- annotate_snps(top_snps_hff)
annotated_hmf <- annotate_snps(top_snps_hmf)

# View results
head(annotated_hff)
head(annotated_hmf)

# Save
write.csv(annotated_hff, "Annotated_HFF_SNPs.csv", row.names = FALSE)
write.csv(annotated_hmf, "Annotated_HMF_SNPs.csv", row.names = FALSE)
```
  git clone https://github.com/MBilgici/Final-Group-Project.git
  cd maize-haploid-gwas-replication

2. Open `code/gwas_hff_hmf_analysis.Rmd` in RStudio
3. Run all chunks or knit the file to HTML
4. Check the `/output` directory for results

### Output
- BLUEs for HFF and HMF traits
- GWAS results from FarmCPU and MLMM
- Manhattan and QQ plots
- Gene annotations of significant SNPs

