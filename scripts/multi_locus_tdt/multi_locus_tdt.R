#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Example script to perform the multi-locus TDT
# using TDT summary statistics and external GWAS summary statistics

# ------------------------------------------------------------
# Input: 
# 1. TDT summary statistics (e.g. tdt_all_off.txt)
# The input file must contain at least the following columns:
#   - SNP      : SNP ID
#   - CHR      : Chromosome
#   - BP       : Base pair position
#   - A1       : Allele 1
#   - A2       : Allele 2
#   - T	       : Transmitted A1 allele count
#   - U	       : Untransmitted A1 allele count
#   - OR       : Odds ratio (allele 1)
#   - P        : TDT p-value

# 2. GWAS summary statistics (e.g. birthweight.txt)
# Summary statistics must be aligned to the allele of interest
# For example, if testing whether birthweight increasing alleles are overtransmitted,
# the EA should be the birthweight increasing allele,
# and BETA should represent the effect size per copy of that allele
# Must contain at least the following columns:
#   - SNP      : SNP ID
#   - EA       : Effect allele
#   - BETA     : GWAS beta (aligned to the EA)
#   - P        : GWAS p-value 

# ------------------------------------------------------------
# Output: 
# Three files for the weighted, unweighted exact and unweighted approximation results 
# (e.g. weighted_approx_multi_loci_tdt_tdt_all_off_birthweight.txt)
# Containing the columns:
#   - SUM              : Observed weighted sum of EA transmissions across all SNPs
#   - EXP              : Expected value of SUM under the null
#   - VAR              : Variance of SUM under the null
#   - Z                : Test statistic
#   - P_normal_approx	 : Multi-locus TDT p-value
#   - n_snps           : Number of SNPs included in the test

# (e.g. unweighted_multi_loci_tdt_tdt_all_off_birthweight.txt)
#   - successes	  : Observed total transmission count across all SNPs
#   - trials	    : Total number of informative transmissions 
#   - pval        : Multi-locus TDT p-value
#   - n_snps      : Number of SNPs included in the test 

# (e.g. unweighted_approx_multi_loci_tdt_tdt_all_off_birthweight.txt)
#   - successes	       : Observed total transmission count across all SNPs
#   - trials	         : Total number of informative transmissions 
#   - SE	             : Standard error
#   - Z	               : Z-score
#   - P_normal_approx  : Multi-locus TDT p-value
#   - n_snps           : Number of SNPs included in the test 

# ------------------------------------------------------------
# Usage:
# for tdt in 'tdt_all_off' 'tdt_female_off' 'tdt_male_off' 'tdt_parent_of_origin'; do
#   for gwas in birthweight.txt 'infertility.txt' ; do
#       Rscript multi_locus_tdt.R ${tdt} ${gwas}
#    done
# done

library(data.table)

args <- commandArgs(trailingOnly=TRUE)

#####################################################
# Read in the datasets
tdt <- fread(paste0(args[1], ".txt"), header=T, select = c("SNP", "CHR", "BP", "A1", "A2", "T", "U", "OR", "P")) 
gwas <- fread(args[2],header=T, select = c("SNP",  "EA", "BETA", "P"))

#####################################################
# Merge GWAS and tdt
merged <- merge(gwas,tdt, by = "SNP")

# Realign the TDT results to reflect the GWAS EA
merged$T_new <- NA
merged$U_new <- NA
merged$T_new <- ifelse(merged$A1 == merged$EA, merged$T, merged$U)
merged$U_new <- ifelse(merged$A1 == merged$EA, merged$U, merged$T)
merged$OR_new <- merged$T_new / merged$U_new

# Subset to columns of interest and rename
merged2 <- merged[,c("SNP", "CHR", "BP", "EA", "BETA", "P.x", "T_new", "U_new", "OR_new", "P.y")]
colnames(merged2) <- c("SNP", "CHR", "BP", "EA", "BETA_GWAS", "P_GWAS", "T", "U", "OR", "P")

#####################################################
# 1. Unweighted multi-locus TDT (exact)
# Count the GWAS EA allele transmissions
# Use the binomial sign test to get the p-value (exact 1 tailed sign test)
# binom.test(x, n, p, alternative="greater")
# x is the number of successes, n is the number of trials, p is the probability of success, one-tailed
successes <- sum(merged2$T) 
trials <- sum(merged2$T) + sum(merged2$U)
pval <- binom.test(x = successes, n = trials, p = 0.5, alternative="greater")$p.value

#####################################################
# 2. Unweighted multi-locus TDT (normal approximation)
# Z = (count - expected_count) / (SE)
# SE = sqrt(np(1-p))
# Right tail
SE = sqrt(trials*0.5^2)
Z = (successes- 0.5*trials) / SE
P_normal_approx = pnorm(q=Z, lower.tail=FALSE)
n_snps = dim(merged2)[1]

# Store the results
unweighted_binomial <- cbind(successes,trials,pval,n_snps)
unweighted_norm_approx <- cbind(successes,trials, SE, Z, P_normal_approx,n_snps)

#####################################################
# 3. Weighted multi-locus TDT (normal approximation)
# Multiply transmission count by the GWAS beta
merged2$weighted_T_count <- merged2$T*merged2$BETA_GWAS
# Sum up the weighted counts
SUM = sum(merged2$weighted_T_count)
# EXP = the expectation, VAR = variance
EXP = 0.5* sum((merged2$T + merged2$U)* merged2$BETA_GWAS)
VAR = 0.25* (sum((merged2$T + merged2$U)*merged2$BETA_GWAS^2))
# Normal approximation
Z = (SUM-EXP)/(sqrt(VAR))
P_normal_approx = pnorm(q=Z, lower.tail=FALSE)
# Store the results
weighted_norm_approx <- cbind(SUM,EXP, VAR, Z, P_normal_approx,n_snps)

# Save the results
write.table(weighted_norm_approx, file=paste0("weighted_approx_multi_loci_tdt_",args[1],"_",args[2]), quote=F, row.names=F, col.names=T, sep="\t")
write.table(unweighted_binomial, file=paste0("unweighted_multi_loci_tdt_",args[1],"_",args[2]), quote=F, row.names=F, col.names=T, sep="\t")
write.table(unweighted_norm_approx, file=paste0("unweighted_approx_multi_loci_tdt_",args[1],"_",args[2]), quote=F, row.names=F, col.names=T, sep="\t")
