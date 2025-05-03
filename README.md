# Final-Group-Project
library(readxl)  
library(dplyr)  
library(CMplot)  
library(ggplot2)  
library(tidyr)  

# Read GWAS Data  
gwas_data <- read_excel("Table 3.XLSX", na = "NA")  

# Create manhattan_data using base R  
cols_needed <- c(2, 5, 8, 13, 14, 15, 16)  
manhattan_data <- gwas_data[, cols_needed]  
colnames(manhattan_data) <- c("SNP", "CHR", "BP", "P_nov_blink", "P_orel_blink", "P_nov_farmcpu", "P_orel_farmcpu")  

# Convert to numeric using base R  
manhattan_data$CHR <- as.numeric(manhattan_data$CHR)  
manhattan_data$BP <- as.numeric(manhattan_data$BP)  
manhattan_data$P_nov_blink <- as.numeric(manhattan_data$P_nov_blink)  
manhattan_data$P_orel_blink <- as.numeric(manhattan_data$P_orel_blink)  
manhattan_data$P_nov_farmcpu <- as.numeric(manhattan_data$P_nov_farmcpu)  
manhattan_data$P_orel_farmcpu <- as.numeric(manhattan_data$P_orel_farmcpu)  

# Remove NA rows  
manhattan_data <- manhattan_data[!is.na(manhattan_data$SNP), ]  

# Create nov_blink using base R  
nov_blink <- data.frame(  
  SNP = manhattan_data$SNP,  
  CHR = manhattan_data$CHR,  
  BP = manhattan_data$BP,  
  P.value = 10^(-as.numeric(manhattan_data$P_nov_blink))  
)  
nov_blink <- nov_blink[!is.na(nov_blink$P.value), ]  

# Create orel_blink  
orel_blink <- data.frame(  
  SNP = manhattan_data$SNP,  
  CHR = manhattan_data$CHR,  
  BP = manhattan_data$BP,  
  P.value = 10^(-as.numeric(manhattan_data$P_orel_blink))  
)  
orel_blink <- orel_blink[!is.na(orel_blink$P.value), ]  

# Create nov_farmcpu  
nov_farmcpu <- data.frame(  
  SNP = manhattan_data$SNP,  
  CHR = manhattan_data$CHR,  
  BP = manhattan_data$BP,  
  P.value = 10^(-as.numeric(manhattan_data$P_nov_farmcpu))  
)  
nov_farmcpu <- nov_farmcpu[!is.na(nov_farmcpu$P.value), ]  

# Create orel_farmcpu  
orel_farmcpu <- data.frame(  
  SNP = manhattan_data$SNP,  
  CHR = manhattan_data$CHR,  
  BP = manhattan_data$BP,  
  P.value = 10^(-as.numeric(manhattan_data$P_orel_farmcpu))  
)  
orel_farmcpu <- orel_farmcpu[!is.na(orel_farmcpu$P.value), ]  

# Calculate significance threshold  
n_markers <- nrow(nov_blink)  
sig_threshold <- 0.05 / n_markers  

# Generate Manhattan plots  
# Novosibirsk BLINK Plot  
CMplot(nov_blink,  
       type = "p",  
       plot.type = c("m","q"),  
       col = c("darkgreen","orange"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "GWAS - Novosibirsk BLINK")  

# Orel BLINK Plot  
CMplot(orel_blink,  
       type = "p",  
       plot.type = c("m","q"),  
       col = c("darkblue","purple"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "GWAS - Orel BLINK")  

# Novosibirsk FarmCPU Plot  
CMplot(nov_farmcpu,  
       type = "p",  
       plot.type = c("m","q"),  
       col = c("brown","magenta"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "GWAS - Novosibirsk FarmCPU")  

# Orel FarmCPU Plot  
CMplot(orel_farmcpu,  
       type = "p",  
       plot.type = c("m","q"),  
       col = c("black","pink"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "GWAS - Orel FarmCPU")  

print("GWAS Manhattan and Q-Q plots generated successfully for all four analyses.")  














library(readxl)  
library(CMplot)  

# ---------------------------  
# 1. Read GWAS Data from Table 3 (FarmCPU results)  
# ---------------------------  
gwas_data <- read_excel("Table 3.XLSX", na = "NA")  

# ---------------------------  
# 2. Prepare common dataframe for Manhattan plots using base R  
# Columns (by index):  
#   Column 2: SNP  
#   Column 3: Trait (DTF or DTM)  
#   Column 5: Chr  
#   Column 8: Pos (BP)  
#   Column 13: BLINK_Novosibirsk (-log10(p))  
#   Column 14: BLINK_Orel (-log10(p))  
#   Column 15: FarmCPU_Novosibirsk (-log10(p))  
#   Column 16: FarmCPU_Orel (-log10(p))  
# ---------------------------  
cols_needed <- c(2, 3, 5, 8, 13, 14, 15, 16)  
manhattan_data <- gwas_data[, cols_needed]  
colnames(manhattan_data) <- c("SNP", "Trait", "CHR", "BP", "P_nov_blink", "P_orel_blink", "P_nov_farmcpu", "P_orel_farmcpu")  

# Convert to numeric where needed and remove NA rows (for SNP, CHR, BP)  
manhattan_data$CHR <- as.numeric(manhattan_data$CHR)  
manhattan_data$BP <- as.numeric(manhattan_data$BP)  
manhattan_data$P_nov_blink <- as.numeric(manhattan_data$P_nov_blink)  
manhattan_data$P_orel_blink <- as.numeric(manhattan_data$P_orel_blink)  
manhattan_data$P_nov_farmcpu <- as.numeric(manhattan_data$P_nov_farmcpu)  
manhattan_data$P_orel_farmcpu <- as.numeric(manhattan_data$P_orel_farmcpu)  
manhattan_data <- manhattan_data[!is.na(manhattan_data$SNP) & !is.na(manhattan_data$CHR) & !is.na(manhattan_data$BP), ]  

# ---------------------------  
# 3. Function to prepare dataset for a given trait and analysis column  
# ---------------------------  
prep_data <- function(df, trait_name, p_column) {  
  subset_df <- df[df$Trait == trait_name, ]  
  # Convert -log10(p) to raw p-value  
  out_df <- data.frame(  
    SNP = subset_df$SNP,  
    CHR = subset_df$CHR,  
    BP  = subset_df$BP,  
    P.value = 10^(-subset_df[[p_column]])  
  )  
  out_df <- out_df[!is.na(out_df$P.value), ]  
  return(out_df)  
}  

# ---------------------------  
# 4. For Trait DTF and DTM, prepare datasets for each analysis:  
# ---------------------------  
# DTF using BLINK (Novosibirsk and Orel) and FarmCPU (Novosibirsk and Orel)  
dtf_nov_blink    <- prep_data(manhattan_data, "DTF", "P_nov_blink")  
dtf_orel_blink   <- prep_data(manhattan_data, "DTF", "P_orel_blink")  
dtf_nov_farmcpu  <- prep_data(manhattan_data, "DTF", "P_nov_farmcpu")  
dtf_orel_farmcpu <- prep_data(manhattan_data, "DTF", "P_orel_farmcpu")  

# DTM using BLINK and FarmCPU  
dtm_nov_blink    <- prep_data(manhattan_data, "DTM", "P_nov_blink")  
dtm_orel_blink   <- prep_data(manhattan_data, "DTM", "P_orel_blink")  
dtm_nov_farmcpu  <- prep_data(manhattan_data, "DTM", "P_nov_farmcpu")  
dtm_orel_farmcpu <- prep_data(manhattan_data, "DTM", "P_orel_farmcpu")  

# ---------------------------  
# 5. Calculate significance threshold based on one dataset sample (e.g., dtf_nov_blink)  
# ---------------------------  
n_markers <- nrow(dtf_nov_blink)  
sig_threshold <- 0.05 / n_markers  

# ---------------------------  
# 6. Generate Manhattan & Q–Q plots using CMplot for each trait and analysis  
# ---------------------------  
# For DTF Trait:  
CMplot(dtf_nov_blink,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("darkgreen", "orange"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTF GWAS - Novosibirsk BLINK")  

CMplot(dtf_orel_blink,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("darkblue", "purple"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTF GWAS - Orel BLINK")  

CMplot(dtf_nov_farmcpu,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("brown", "magenta"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTF GWAS - Novosibirsk FarmCPU")  

CMplot(dtf_orel_farmcpu,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("black", "pink"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTF GWAS - Orel FarmCPU")  

# For DTM Trait:  
CMplot(dtm_nov_blink,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("darkgreen", "orange"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTM GWAS - Novosibirsk BLINK")  

CMplot(dtm_orel_blink,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("darkblue", "purple"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTM GWAS - Orel BLINK")  

CMplot(dtm_nov_farmcpu,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("brown", "magenta"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTM GWAS - Novosibirsk FarmCPU")  

CMplot(dtm_orel_farmcpu,  
       type = "p",  
       plot.type = c("m", "q"),  
       col = c("black", "pink"),  
       threshold = sig_threshold,  
       threshold.lty = 1,  
       threshold.col = "red",  
       threshold.lwd = 1,  
       file.output = FALSE,  
       main = "DTM GWAS - Orel FarmCPU")  

cat('GWAS Manhattan and Q-Q plots generated successfully for traits DTF and DTM using both models.\n')  

