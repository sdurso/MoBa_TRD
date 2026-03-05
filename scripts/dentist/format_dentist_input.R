#!/usr/bin/env Rscript
# Example script to format TDT summary statistics for DENTIST

# Input: PLINK TDT output file with additional columns e.g. tdt_all_off.txt
# The input file must contain at least the following columns (from PLINK TDT output):
#   - SNP      : SNP ID
#   - CHR      : Chromosome
#   - BP       : Base pair position
#   - A1       : Allele 1
#   - A2       : Allele 2
#   - T	       : Transmitted A1 allele count
#   - U	       : Untransmitted A1 allele count
#   - OR       : Odds ratio (allele 1)
#   - CHISQ	   : TDT chi-square statistic
#   - P        : TDT p-value
# Additional columns: 
#   - Total    : Total number of informative trios (T + U)
#   - MAF      : Frequency of allele 1

# Output: dentist_*_chr*.txt
# Reformatted summary statistics for DENTIST containing the columns:
#   - SNP      : SNP ID
#   - A1       : Allele 1
#   - A2       : Allele 2
#   - freq     : Frequency of allele 1
#   - beta     : Z-score
#   - se	     : Standard error of the Z-score (1)
#   - p        : TDT p-value
#   - N        : Sample size

# Usage:
# for i in 'tdt_all_off' 'tdt_female_off' 'tdt_male_off' 'tdt_parent_of_origin'; do
#     Rscript format_dentist_input.R ${i}
# done

#####################################################
library(data.table)

# Get the TDT base name from command line
args <- commandArgs(trailingOnly=TRUE)
tdt_base <- args[1]

#####################################################
# Read TDT summary statistics
data <- fread(paste0(tdt_base, ".txt"))

# Compute beta (Z-scores)
# Make beta negative when OR < 1, otherwise leave it positive
data$beta <- ifelse(data$OR < 1, -sqrt(data$CHISQ), sqrt(data$CHISQ))
# Set se as 1
data$se <- 1

# Rename columns
# DENTIST requires these specific columns
data$N <- data$Total
data$freq <- data$MAF
data$p <- data$P

# Split by chromosome
chromosomes <- unique(data$CHR)

# Save per chromosome
for (chr in chromosomes) {
  subset_data <- data[data$CHR == chr, ]
  tmp <- subset_data[order(subset_data$BP),]
  tmp2 <- tmp[,c("SNP", "A1", "A2", "freq", "beta", "se", "p", "N")]
  filename <- paste0("dentist_", tdt_base, "_chr", chr, ".txt")
  write.table(tmp2, file=filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ")
}
