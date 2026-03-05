#!/usr/bin/env Rscript
# Rscript that will simulates trio genotypes under transmission ratio distortion (TRD)
# and selection to evaluate how selection affects Hardy–Weinberg equilibrium (HWE)
# and transmission disequilibrium test (TDT) test statistics


# Usage:
#   Rscript simulate_selection_trd_tdt_hwe.R 0.3 0.05
# Where arguments are:
#  - MAF : minor allele frequency)
#  - s   : selection coefficient)

# Output:
#   e.g. sim_results_hwe_tdt_p_0.3_s_0.05.txt

##############################################################################################################################################
library(data.table)
library(gaston)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

# First arg is p
# Second arg is s
# e.g.
#args <- c(0.5, 0.1)
Sys.time()
set.seed(27)

##############################################################################################################################################
# Simulation parameters 
N <- 43841 # Number of families
Nreps <- 1000 # Number of replications
####################################################################################################################################################
# Set up vectors and matrices to store variables:
# For each SNP:
MU <- vector(mode = "numeric", length = N) # Maternal untransmitted allele across all individuals for one SNP
MT <- vector(mode = "numeric", length = N) # Maternal transmitted allele for one SNP
PU <- vector(mode = "numeric", length = N) # Paternal untransmitted allele for one SNP
PT <- vector(mode = "numeric", length = N) # Paternal transmitted allele for one SNP
M_geno <- vector(mode = "numeric", length = N) # Maternal genotypes
P_geno <- vector(mode = "numeric", length = N) # Paternal genotypes
O_geno <- vector(mode = "numeric", length = N) # Offspring genotypes

#######################################################################################################################################
# Function for generating parental genotypes and offspring genotypes
# Later I select upon offspring genotype
# Requires MAF and amount of transmission distortion as the input (default TRD is zero)
transmissions <- function(p_input, td_input_mat = 0, td_input_pat = 0, s_input, h_input){
  off_genotypes <- c()
  mat_genotypes <- c()
  pat_genotypes <- c()
  mat_transmissions <- c()
  pat_transmissions <- c()
  mat_zygosity <- c()
  pat_zygosity <- c()
  off_zygosity <- c()
  
  for (j in 1:Nreps){
    # Set the allele frequencies and amount of transmission distortion
    p_snp <- p_input
    q_snp <- 1-p_snp
    TD_snp_mat <- td_input_mat
    TD_snp_pat <- td_input_pat
    
    # First sample the maternal and paternal genotypes (based on HWE) for each mother and father
    M_geno <- sample(x = c(0,1,2), replace = T, size = N, prob=c(p_snp*p_snp, 2*p_snp*q_snp, q_snp*q_snp))
    P_geno <- sample(x = c(0,1,2), replace = T, size = N, prob=c(p_snp*p_snp, 2*p_snp*q_snp, q_snp*q_snp))
    
    # For each family
    # Simulate the transmitted and non-transmitted allele
    # When there is TD, the allele transmitted to offspring from heterozygous parent is not equally likely
    # When there is selection, the allele transmitted to offspring from parent is not equally likely
    
    for(i in 1:N){ 
      if(M_geno[i] == 0){ # If the mother is homozygous for the SNP (aa)
        MT[i] <- 0 # She must transmit allele a to offspring
        MU[i] <- 0 # The non-transmitted allele is also a
      } else if (M_geno[i] == 2){ # Same applies for a mother homozygous for the SNP (AA)
        MT[i] <- 1
        MU[i] <- 1
      } else if (M_geno[i] == 1) { # If the mother is heterozygous, then either allele could be transmitted
        # Make it more likely that allele a is transmitted 
        # Do this by altering the probability from 0.5 to something greater!
        MT[i] <- sample(x = c(0,1), replace = T, size = 1, prob=c(0.5+TD_snp_mat,0.5-TD_snp_mat))
        MU[i] <- M_geno[i] - MT[i]
      }
      
      # Same for fathers
      if(P_geno[i] == 0){
        PT[i] <- 0
        PU[i] <- 0
      }else if (P_geno[i] == 2){
        PT[i] <- 1
        PU[i] <- 1
      } else if (P_geno[i] == 1) {
        PT[i] <- sample(x = c(0,1), replace = T, size = 1, prob=c(0.5+TD_snp_pat,0.5-TD_snp_pat))
        # If i wanted the paternal effect to be in the opposite direction and cancel out maternal effect
        # PT[i] <- sample(x = c(0,1), replace = T, size = 1, prob=c(0.5-TD,0.5+TD)) 
        PU[i] <- P_geno[i] - PT[i]
      }
    }
    # Set the offspring genotype
    O_geno <- MT + PT
    
    # Simulate survival
    # s is the selection coefficient (proportion of fitness lost due to a particular genotype)
    # h is the measure of dominance (h = 0.5 means that it is additive)
    # The homozygote (0) is not selected against
    # The heterozygote (1) is selected against
    # The homozygote (2) is selected against
    
    fitness_aa = 1
    fitness_Aa = 1-h_input*s_input
    fitness_AA = 1-s_input
    
    # Simulate values between 0 and 1 for each offspring
    temp <- runif(n = N, min = 0, max = 1)
    
    # Now perform the selection
    # Selection  is acting against the A allele 
    # If the offspring has genotype AA (2), and their random number is greater than fitness_AA, then exclude them
    # if s is 0.9, then 90% of the AA offspring should become NA
    # Also NA the parental transmissions
    
    for(i in 1:N){ 
      if (O_geno[i] == 2 & temp[i] > fitness_AA){ 
        O_geno[i] <- NA
        MT[i] <- NA
        PT[i] <- NA
      } else if (O_geno[i] == 1 & temp[i] > fitness_Aa) { 
        O_geno[i] <- NA
        MT[i] <- NA
        PT[i] <- NA
      }
    }

    # Bind the vectors into a matrix
    off_genotypes <- cbind(off_genotypes, O_geno)
    mat_genotypes <- cbind(mat_genotypes, M_geno)
    pat_genotypes <- cbind(pat_genotypes, P_geno)
    
    # Bind the vectors into a matrix
    mat_transmissions <- cbind(mat_transmissions, MT)
    pat_transmissions <- cbind(pat_transmissions, PT)
    mat_zygosity <- cbind(mat_zygosity, M_geno)
    pat_zygosity <- cbind(pat_zygosity, P_geno)
    off_zygosity <- cbind(off_zygosity, O_geno)
  }
  
  # Name the replications
  rep_name <- paste0(c("Rep_"), seq(1:Nreps))
  colnames(off_genotypes) <- rep_name
  
  # Add a column for offspring ID and SEX (2 = female, 1 = male)
  ID <- paste0(c("Offspring_"), seq(1:N))
  off_genotypes <- cbind(ID, off_genotypes)
  
  # Return these results
  out <- list("off_genotypes" = off_genotypes, "mat_transmissions" = mat_transmissions, "pat_transmissions" = pat_transmissions, "mat_zygosity" = mat_zygosity, "pat_zygosity"= pat_zygosity, "off_zygosity"= off_zygosity )
  
}

#######################################################################################################################################
# Function to perform the TDT
tdt <- function(maternal_transmissions, maternal_zygosity, paternal_transmissions, paternal_zygosity){
  # Setup TDT results tables
  tdt_both <- matrix(data = NA, nrow = Nreps, ncol = 5)
  colnames(tdt_both) <- c("Total", "A1_count", "A2_count", "X_sq", "P")
  rep_name <- paste0(seq(1:Nreps))
  rownames(tdt_both) <- rep_name
  
  tdt_maternal <- matrix(data = NA, nrow = Nreps, ncol = 5)
  colnames(tdt_maternal) <- c("M_Total", "M_A1_count", "M_A2_count", "M_X_sq", "M_P")
  rownames(tdt_maternal) <- rep_name
  
  tdt_paternal <- matrix(data = NA, nrow = Nreps, ncol = 5)
  colnames(tdt_paternal) <- c("F_Total", "F_A1_count", "F_A2_count", "F_X_sq", "F_P")
  rownames(tdt_paternal) <- rep_name
  
  colnames(maternal_zygosity) <- rep_name
  colnames(paternal_zygosity) <- rep_name
  colnames(maternal_transmissions) <- rep_name
  colnames(paternal_transmissions) <- rep_name
  
  for(i in seq(length(rep_name))){  # Iterate through each SNP
    rep <- rep_name[i]
    
    # Extract only offspring whose parent/parents are heterozygous
    # Then count the transmissions of the minor and major allele
    # Order of SNPs and individuals are conserved across data
    
    # Subset the maternal_transmission to only contain offspring whose mothers were heterozygous for SNP i
    hets_m <- subset(maternal_transmissions, maternal_zygosity[,paste0(rep)] == 1)
    sum_A1_m <- length(which(hets_m[,paste0(rep)]==0))
    sum_A2_m <- length(which(hets_m[,paste0(rep)]==1))
    total_m <- nrow(hets_m)
    
    # For paternal transmissions
    hets_p <- subset(paternal_transmissions, paternal_zygosity[,paste0(rep)] == 1)
    sum_A1_p <- length(which(hets_p[,paste0(rep)]==0)) 
    sum_A2_p <- length(which(hets_p[,paste0(rep)]==1)) 
    total_p <- nrow(hets_p)
    
    # TDT
    sum_A1_both <- sum_A1_m + sum_A1_p 
    sum_A2_both <- sum_A2_m + sum_A2_p
    total_both <- total_m + total_p
    chi_both <- (sum_A1_both - sum_A2_both)^2 / (sum_A1_both + sum_A2_both)
    p_both <- pchisq(chi_both, df=1, lower.tail=FALSE)
    tdt_both[i,"A1_count"] <- sum_A1_both
    tdt_both[i,"A2_count"] <- sum_A2_both
    tdt_both[i,"Total"] <- total_both
    tdt_both[i,"X_sq"] <- chi_both
    tdt_both[i,"P"] <- p_both 
  }
  # Return the results
  out <- list("tdt_both" = tdt_both, "tdt_maternal" = tdt_maternal, "tdt_paternal" = tdt_paternal)
  return(out)
}

###########################################################################################
# Perform the HWE Chisq test

# Set up the HWE chisq test function
# requires offspring genotypes and amount of transmission distortion as the input
hwe <- function(off_geno_input, td_input_mat = 0 , td_input_pat = 0 , s_input, h_input){
  # Setup HWE test result table
  rep_name <- paste0(c("Rep_"), seq(1:Nreps))
  HWE <- matrix(data = NA, nrow = length(rep_name), ncol = 10)
  colnames(HWE) <- c("a_freq", "A_freq", "aa_obs","Aa_obs","AA_obs","aa_exp","Aa_exp","AA_exp", "Chi", "P")
  rownames(HWE) <- rep_name
  
  for(i in seq(length(rep_name))){  # Iterate through each replication
    rep <- rep_name[i]
    
    # Observed genotype counts
    aa_obs <- length(which(off_geno_input[,paste0(rep)] == 0))
    Aa_obs <- length(which(off_geno_input[,paste0(rep)] == 1))
    AA_obs <- length(which(off_geno_input[,paste0(rep)] == 2))
    
    # Observed genotype frequencies
    # N = number of offspring
    aa_freq <- aa_obs/N
    Aa_freq <- Aa_obs/N
    AA_freq <- AA_obs/N
    
    # Observed allele frequencies
    a_freq <- (2*aa_obs + Aa_obs)/(2*N)
    A_freq <- (2*AA_obs + Aa_obs)/(2*N)
    
    # Expected genotype counts
    aa_exp <- N*(a_freq^2)  # N x p^2
    Aa_exp <- N*(2*a_freq*A_freq) # N x 2pq
    AA_exp <- N*(A_freq^2) # N x q^2
    
    # Perform the chisq test
    chi <- (((aa_obs - aa_exp)^2)/aa_exp) + (((Aa_obs - Aa_exp)^2)/Aa_exp) + (((AA_obs - AA_exp)^2)/AA_exp)
    # Under the null the chisq test statistic has an approximate chisq distr with 1df
    # Mean Chi is the number of degrees of freedom
    p <- pchisq(chi, df=1, lower.tail=FALSE)
    
    # Sub the values into the results table
    HWE[i,"a_freq"] <- a_freq
    HWE[i,"A_freq"] <- A_freq
    HWE[i,"aa_obs"] <- aa_obs
    HWE[i,"Aa_obs"] <- Aa_obs
    HWE[i,"AA_obs"] <- AA_obs
    HWE[i,"aa_exp"] <- aa_exp
    HWE[i,"Aa_exp"] <- Aa_exp
    HWE[i,"AA_exp"] <- AA_exp
    HWE[i,"Chi"] <- chi
    HWE[i,"P"] <- p
  }
  # New column containing the amount of TD
  t_m <- rep(td_input_mat, times = Nreps)
  t_p <- rep(td_input_pat, times = Nreps)
  s <- rep(s_input, times = Nreps)
  h <- rep(h_input, times = Nreps)
  HWE <- cbind(HWE, t_m, t_p, s, h)
  # Return the results
  return(HWE)
}

###########################################################################################
# Simulate various scenarios - all with no TRD
transmissions_a <- transmissions(p_input = as.numeric(args[1]),s_input= as.numeric(args[2]), h_input =0.5)
tdt_a <- tdt(maternal_transmissions = transmissions_a$mat_transmissions, maternal_zygosity= transmissions_a$mat_zygosity, 
             paternal_transmissions = transmissions_a$pat_transmissions, paternal_zygosity = transmissions_a$pat_zygosity)
HWE_a <- hwe(off_geno_input = transmissions_a$off_genotypes,s_input= as.numeric(args[2]), h_input =0.5)
print("a complete")

a <- as.data.table(cbind(HWE_a[,c("P")], tdt_a$tdt_both[,c("P")], rep(as.numeric(args[2]), times=Nreps), rep(as.numeric(args[1]), times=Nreps)))
colnames(a) <- c("HWE_P", "TDT_P", "S", "MAF")
a$HWE_P_log <- -log10(a$HWE_P)
a$TDT_P_log <- -log10(a$TDT_P)
a$REP <- rownames(HWE_a)
write.table(a, file = paste0("sim_results_hwe_tdt_p_", args[1],"_s_", args[2],".txt"), quote=F, row.names=F, col.names=T)
print("table written")
