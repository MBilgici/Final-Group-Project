# Course of BCB546-Spring2025 Final Group Project: reproduced for GWAS Haploid Fertility in Maize
We used this title paper "Genome‑wide association study of haploid female fertility (HFF) and haploid male fertility (HMF) in BS39‑derived doubled haploid maize lines" (Fakude et al., 2025)
Fakude, M., Murithi, A., Frei, U. K., Scott, P. M., & Lübberstedt, T. (2025). Genome-wide association study of haploid female fertility (HFF) and haploid male fertility (HMF) in BS39-derived doubled haploid maize lines. Theoretical and Applied Genetics, 138(1), 1-14.
## Overview
This repository contains all the materials related to the reproducing of the following study:

Our target is to reproduce and validate the GWAS analyses presented in the orginal paper using the original data and methods, and to document all findings in a reproducible format.

## Repository Structure
1-README.md
2-Fakude et al., 2025.md (Summary and analysis of the original (Fakude et al., 2025) paper)
3-R code (R scripts and RMarkdown for data processing and analysis)
4-data data sources
5-Our results output (Generated plots, summary tables, and annotations) GWAS

## When we reproduced the original paper and follewed the below steps
Install & load everything you need
Read your HFF & HMF CSVs, compute binomial‐GLM BLUES for each
Merge those into one phenotype file
Patch & filter your HapMap to biallelic SNPs
Run rMVP (FarmCPU & MLM) on both traits at once
Produce Manhattan & QQ plots for each trait
Extract genome‐wide significant SNPs and annotate HMF hits via Ensembl Plants
### Requirements
R and R studio Software
Data from Original Paper (you can find our github also)

### Instructions
1. Clone this repository:
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


## HFF Data Processing and BLUEs
## GWAS Analysis with RVRW

setwd("C:/Users/bilgici/Desktop/New Folder/2025")

# 1) install & load packages

pkgs <- c("lme4","ggResidpanel","dplyr","ggplot2","emmeans","tidyverse",
          "car","qqman","readxl","rMVP","data.table","bigmemory","arm")
for(p in pkgs){
  if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
  library(p, character.only=TRUE)
}

# 2) Read & process HMF data, fit GLM, get emmeans (BLUES) ------------------
da <- read.csv("HFF.csv")
da$PI   <- as.factor(da$PI)
da$Rep  <- as.factor(da$Rep)
da$HFF  <- as.numeric(da$HFF)

# Check distribution if you like
hist(da$HFF, col="cyan", xlab="Haploid male fertility")

# Aggregate Fertile/Sterile counts by PI × Rep
da.agg <- da %>%
  group_by(PI,Rep) %>%
  summarize(
    Fertile = sum(Fertile),
    Sterile = sum(Sterile),
    .groups="drop"
  ) %>%
  mutate(
    HFF = Fertile/(Fertile+Sterile),
    HFF.mod = ifelse(HFF==0,
                     min(HFF[HFF>0])/2,
                     HFF)
  )

# Fit binomial GLM
glm_GWAS <- glm(
  cbind(Fertile,Sterile) ~ Rep + PI,
  data = da.agg,
  family = binomial
)
resid_panel(glm_GWAS)

# Extract emmeans on the PI effect and back‐transform to %  
emm <- emmeans(glm_GWAS, "PI")
glm_eblues_GWAS <- data.frame(emm)
glm_eblues_GWAS$HFF <- arm::invlogit(glm_eblues_GWAS$emmean)*100

pheno_rMVP <- glm_eblues_GWAS %>%
  dplyr::select(Taxa = PI, HFF)
pheno_rMVP <- glm_eblues_GWAS[, c("PI", "HFF")]
names(pheno_rMVP) <- c("Taxa", "HFF")

# write out
write.csv(pheno_rMVP, "HFF_BLUES.csv", row.names=FALSE, quote=FALSE)

# 3) Patch and filter your genotype HapMap as before ------------------------
lines <- readLines("genotype_data.txt")
lines[1] <- sub("^(rs)(\t)", "rs#\\2", lines[1])
writeLines(lines, "genotype_data_fixed.hmp.txt")

hap <- data.table::fread("genotype_data_fixed.hmp.txt", data.table=FALSE)
biallelic <- hap[
  sapply(gregexpr("/", hap$alleles), function(x) sum(x>0)) == 1 &
    sapply(strsplit(hap$alleles, "/"), function(a) all(nchar(a)==1)),
]
cat(sprintf(
  "Kept %d biallelic SNPs (removed %d multi-allelic)\n",
  nrow(biallelic), nrow(hap)-nrow(biallelic)
))
write.table(
  biallelic,
  file = "genotype_data_biallelic.hmp.txt",
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names=TRUE
)

# 4) Make results folder ---------------------------------------------------
if(!dir.exists("rMVP_results")) dir.create("rMVP_results")

# 5) Prepare all data for GWAS with rMVP -------------------------------
MVP.Data(
  fileHMP    = "genotype_data_biallelic.hmp.txt",
  filePhe    = "HFF_BLUES.csv",
  sep.phe    = ",",
  fileKin    = TRUE,
  filePC     = TRUE,
  SNP.impute = "Major",
  out        = "GWAS"
)

# 6) Attach matrices -------------------------------------------------------
map  <- read.table("GWAS.geno.map", header=TRUE, sep="\t", stringsAsFactors=FALSE)
phe  <- read.table("GWAS.phe",     header=TRUE, sep="\t", stringsAsFactors=FALSE)
geno <- attach.big.matrix("GWAS.geno.desc")
K    <- attach.big.matrix("GWAS.kin.desc")
PC   <- attach.big.matrix("GWAS.pc.desc")

# 7) Run FarmCPU & MLM ------------------------------------------------------
res <- MVP(
  phe         = phe,
  geno        = geno,
  map         = map,
  K           = K,
  CV.GLM      = PC,
  CV.MLM      = PC,
  CV.FarmCPU  = PC,
  method      = c("FarmCPU","MLM"),
  nPC.GLM     = 0,
  nPC.MLM     = 0,
  nPC.FarmCPU = 0,
  file.output = TRUE,
  outpath     = "rMVP_results"
)

# Script complete: check rMVP_results/ for HFF results.


## HMF Data Processing and BLUEs


2) Read & process HMF data, fit GLM, get emmeans (BLUES) ------------------
da <- read.csv("HMF.csv")
da$PI   <- as.factor(da$PI)
da$Rep  <- as.factor(da$Rep)
da$HMF  <- as.numeric(da$HMF)

# Check distribution if you like
hist(da$HMF, col="cyan", xlab="Haploid male fertility")

# Aggregate Fertile/Sterile counts by PI × Rep
da.agg <- da %>%
  group_by(PI,Rep) %>%
  summarize(
    Fertile = sum(Fertile),
    Sterile = sum(Sterile),
    .groups="drop"
  ) %>%
  mutate(
    HMF = Fertile/(Fertile+Sterile),
    HMF.mod = ifelse(HFF==0,
                     min(HMF[HMF>0])/2,
                     HMF)
  )

# Fit binomial GLM
glm_GWAS <- glm(
  cbind(Fertile,Sterile) ~ Rep + PI,
  data = da.agg,
  family = binomial
)
resid_panel(glm_GWAS)

# Extract emmeans on the PI effect and back‐transform to %  
emm <- emmeans(glm_GWAS, "PI")
glm_eblues_GWAS <- data.frame(emm)
glm_eblues_GWAS$HMF <- arm::invlogit(glm_eblues_GWAS$emmean)*100

pheno_rMVP <- glm_eblues_GWAS %>%
  dplyr::select(Taxa = PI, HMF)
pheno_rMVP <- glm_eblues_GWAS[, c("PI", "HMF")]
names(pheno_rMVP) <- c("Taxa", "HMF")



# write out
write.csv(pheno_rMVP, "HMF_BLUES.csv", row.names=FALSE, quote=FALSE)

# 3) Patch and filter your genotype HapMap as before ------------------------
lines <- readLines("genotype_data.txt")
lines[1] <- sub("^(rs)(\t)", "rs#\\2", lines[1])
writeLines(lines, "genotype_data_fixed.hmp.txt")

hap <- data.table::fread("genotype_data_fixed.hmp.txt", data.table=FALSE)
biallelic <- hap[
  sapply(gregexpr("/", hap$alleles), function(x) sum(x>0)) == 1 &
    sapply(strsplit(hap$alleles, "/"), function(a) all(nchar(a)==1)),
]
cat(sprintf(
  "Kept %d biallelic SNPs (removed %d multi-allelic)\n",
  nrow(biallelic), nrow(hap)-nrow(biallelic)
))
write.table(
  biallelic,
  file = "genotype_data_biallelic.hmp.txt",
  sep = "\t", quote = FALSE,
  row.names = FALSE, col.names=TRUE
)

# 4) Make results folder ---------------------------------------------------
if(!dir.exists("rMVP_results")) dir.create("rMVP_results")

# 5) Prepare all data for GWAS with rMVP -------------------------------
MVP.Data(
  fileHMP    = "genotype_data_biallelic.hmp.txt",
  filePhe    = "HMF_BLUES.csv",
  sep.phe    = ",",
  fileKin    = TRUE,
  filePC     = TRUE,
  SNP.impute = "Major",
  out        = "GWAS"
)

# 6) Attach matrices -------------------------------------------------------
map  <- read.table("GWAS.geno.map", header=TRUE, sep="\t", stringsAsFactors=FALSE)
phe  <- read.table("GWAS.phe",     header=TRUE, sep="\t", stringsAsFactors=FALSE)
geno <- attach.big.matrix("GWAS.geno.desc")
K    <- attach.big.matrix("GWAS.kin.desc")
PC   <- attach.big.matrix("GWAS.pc.desc")

# 7) Run FarmCPU & MLM ------------------------------------------------------
res <- MVP(
  phe         = phe,
  geno        = geno,
  map         = map,
  K           = K,
  CV.GLM      = PC,
  CV.MLM      = PC,
  CV.FarmCPU  = PC,
  method      = c("FarmCPU","MLM"),
  nPC.GLM     = 0,
  nPC.MLM     = 0,
  nPC.FarmCPU = 0,
  file.output = TRUE,
  outpath     = "rMVP_results"
)

# Heritibility anaysis 
library(lme4)
library(ggResidpanel)
library(dplyr)
library(ggplot2)
library(emmeans)
library(arm)
library(tidyverse)
library(MASS)
library(bestNormalize)


#Read in HFF data

data1 <- read.csv("./HFF.csv")
data1$PI = as.factor(data1$PI)
data1$Rep = as.factor(data1$Rep)
data1$HFF = as.numeric(data1$HFF)


hist(data1$HFF)
hist(data1$HFF, col = "Cyan", xlab = "Haploid Female Fertility", main = "Trait distribution")
Boxplot(data1$HFF, col = "Cyan", ylab = "Haploid Female Fertility", main = "Haploid Female Fertility")


# Check for assumptions: P-Value less than 0.05 violates assumption. Can be seen with histogram
shapiro.test(data1$HFF) 


# the code below will find a transformation that works for your data
bestNormalize::bestNormalize(data1$HFF) # read its output
tmp <- bestNormalize::bestNormalize(data1$HFF)
shapiro.test(tmp$chosen_transform$x.t) # best transformation
data1$transformed_data <- tmp$chosen_transform$x.t # adding transformed data
hist(data1$transformed_data) # better?
hist(data1$transformed_data, col = "Cyan", xlab = "Haploid Female Fertility", main = "Trait distribution")

# ANOVA model
str(data1)
table(data1$PI, data1$Rep)
model <- lm(transformed_data~Rep+PI, data = data1)
#model <- lmer(HFF ~ Genotype + (1 | Rep), data = data1) 

ggResidpanel::resid_panel(model)
anova(model) 

#Calculate Entry-mean heritability from ANOVA output
#PI_variance = Genetic variance
#Genetic variance= (PI MSquare - Residual MSquare)/2
heritability <- PI_variance / (PI_variance + (residual_variance / n))


#For HMF Entry-mean heritability


library(lme4)
#Read in HFF data
da <- read.csv("./HMF.csv")

# Install and load necessary package
install.packages("lme4")
library(lme4)

# Fit a binomial GLMM
glmm_model <- glmer(cbind(Fertile, Sterile) ~ (1 | PI) + (1 | Rep), 
                    family = binomial, data = da.agg)

# Print the summary of the model
summary(glmm_model)
resid_panel(glmm_model)

# Extract variance components
varcomp <- as.data.frame(VarCorr(glmm_model))

# Variance for PI
pi_variance <- varcomp[varcomp$grp == "PI", "vcov"]

# Variance for Rep
rep_variance <- varcomp[varcomp$grp == "Rep", "vcov"]

# Calculate residual variance
deviance_residuals <- residuals(glmm_model, type = "deviance")
residual_variance <- var(deviance_residuals)

# Print variance components
cat("Variance for PI:", pi_variance, "\n")
cat("Variance for Rep:", rep_variance, "\n")
cat("Residual Variance:", residual_variance, "\n")


# Variance components 
pi_variance <- 2.172924
residual_variance <- 4.572366
rep_variance <- 0.03924474

# Number of environments 
n <- 2

# Calculate entry-mean heritability
heritability <- pi_variance / (pi_variance + (residual_variance / n))

# Print the heritability
cat("Entry-Mean Heritability (H²):", heritability, "\n")



 ##
 git clone https://github.com/MBilgici/Final-Group-Project.git
  

2. Open `code/Group Project 2025.Rmd` in RStudio
3. Run all chunks or knit the file to HTML
4. Check the `/output` directory for results

### Output
- BLUEs for HFF and HMF traits
- GWAS results from FarmCPU and MLMM
- Manhattan and QQ plots
- Gene annotations of significant SNPs

## Results for graphs and images
