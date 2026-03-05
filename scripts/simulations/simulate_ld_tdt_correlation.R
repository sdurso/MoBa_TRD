#!/usr/bin/env Rscript
# Rscript that will simulate data
# Test whether the correlation between TDT Z-scores matches the genetic LD under various conditions

# Output: 
# 1. simulation_results.txt
# 2. tdt_ld_correlation.png

library(dplyr)
library(ggplot2)

#####################################################
# 1. Function that simulate two SNP haplotypes in LD given MAF and r2

simulate_ld_haplotypes <- function(
  n_ind, # number of individuals
  maf1 = 0.3, # maf1, maf2: minor allele frequencies for loci 1 and 2
  maf2 = 0.3,
  r2   = 0.5, # linkage disequilibrium measure (0 = independent, 1 = perfect correlation)
  seed = NULL 
) {
  if (!is.null(seed)) set.seed(seed)

  pA <- maf1
  pB <- maf2
  
  # Dmax is the maximum possible linkage disequilibrium coefficient given the allele frequencies
  Dmax <- min(pA * (1 - pB), (1 - pA) * pB)
  
  # Computes the covariance between alleles
  D <- sqrt(r2 * pA * (1 - pA) * pB * (1 - pB))

  if (D > Dmax) stop("Requested r2 not achievable.")
  
  # Compute haplotype frequencies
  # Correlated alleles = +D
  hap_freqs <- c(
    AB = pA * pB + D,
    Ab = pA * (1 - pB) - D,
    aB = (1 - pA) * pB - D,
    ab = (1 - pA) * (1 - pB) + D
  )
  
  # Sample haplotypes (two per individual) using the computed frequencies
  haps <- sample(
    names(hap_freqs),
    size = 2 * n_ind,
    replace = TRUE,
    prob = hap_freqs
  )
  
  # Returns a matrix where each row contains one individual's two haplotypes
  matrix(haps, nrow = n_ind, ncol = 2, byrow = TRUE)
}

#####################################################
# 2. Function that simulates Mendelian transmission 
# of one haplotype from a parent to offspring,
# with optional transmission distortion

# input = a parent's 2 haplotypes
# p_transmit = probability of transmitting 'A' allele (0.5 = Mendelian)
transmit_haplotype <- function(haps, p_transmit = 0.5) {
  # Extract first allele from each haplotype
  allele1 <- substr(haps, 1, 1)
  
  # If homozygous at locus 1, randomly pick one to transmit
  if (allele1[1] == allele1[2]) {
    return(sample(haps, 1))
  }
  
  # If heterozygous
  # Find which position has the A allele
  idx_A <- which(allele1 == "A")
  # Generate a random number between 0 and 1
  # If < p_transmit, transmit A
  # Otherwise transmit a
  if (runif(1) < p_transmit) haps[idx_A] else haps[-idx_A]
}

#####################################################
# 3.Function to simulate parent-offspring trio haplotypes 
# with optional transmission distortion

simulate_trios <- function(
  n_trios, # number of trios
  maf1 = 0.3, # maf1, maf2: minor allele frequencies for loci 1 and 2
  maf2 = 0.3,
  r2   = 0.5, # linkage disequilibrium between the two loci
  p_transmit = 0.5, # probability of transmitting A allele
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Generate parental haplotypes under LD
  fathers <- simulate_ld_haplotypes(n_trios, maf1, maf2, r2)
  mothers <- simulate_ld_haplotypes(n_trios, maf1, maf2, r2)
  
  # Simulate transmission from each parent to child with optional TRD
  child <- matrix(NA, n_trios, 2)
  for (i in 1:n_trios) {
    child[i, 1] <- transmit_haplotype(fathers[i, ], p_transmit)
    child[i, 2] <- transmit_haplotype(mothers[i, ], p_transmit)
  }

  # Label columns and return as list
  colnames(fathers) <- colnames(mothers) <- colnames(child) <- c("hap1", "hap2")
  list(father = fathers, mother = mothers, child = child)
}

#####################################################
# 4. Function that identifies which allele was transmitted 
# from parent to child at a specific SNP


# parent_haps: parent's two haplotypes (e.g. 'AB, 'ab'))
# child_haps: child's two haplotypes (one from each parent)
# snp: which locus to check (1 or 2)
get_transmitted_allele <- function(parent_haps, child_haps, snp) {
  # Loop over parent haplotypes
  for (h in parent_haps) {
    # If the child has this haplotype, consider it transmitted
    if (sum(child_haps == h) > 0) {
      return(substr(h, snp, snp))
    }
  }
  stop("No transmitted haplotype found from this parent!")
}

#####################################################
#5. Function that counts transmitted vs untransmitted risk alleles 
# from heterozygous parents at a single SNP

# trios: output from simulate_trios() containing father, mother, child haplotypes
# snp: which locus to test (1 or 2)
# risk_allele: which allele to count as 'risk'
compute_tdt <- function(trios, snp = 1, risk_allele = "A") {
  # T = transmitted count
  # U = untransmitted count
  # Both start at zero
  T <- 0
  U <- 0

  # Loop through each trio
  for (i in 1:nrow(trios$child)) {
    for (parent in c("father", "mother")) {
      # Extract parental and child haplotypes
      parent_haps <- trios[[parent]][i, ]
      child_haps <- trios$child[i, ]
      
      # Extract alleles at the SNP of interest from parent's two haplotypes
      parent_alleles <- substr(parent_haps, snp, snp)
      
      # Only informative if parent heterozygous
      if (parent_alleles[1] != parent_alleles[2]) {
        # Determine which allele was transmitted to the child
        transmitted <- get_transmitted_allele(parent_haps, child_haps, snp)
        # Increase the count according to which allele was transmitted
        if (transmitted == risk_allele) T <- T + 1 else U <- U + 1
      }
    }
  }
  # McNemar's test
  z <- (T - U) / sqrt(T + U)
  list(T = T, U = U, Z = z, chi_sq = z^2)
}

#####################################################
# 6. Function that conducts the TDT across multiple SNPs

# trios = haplotype data for the trio
# risk_alleles = vector specifying which allele is the risk one at each SNP
compute_tdt_multi <- function(trios, risk_alleles) {
  # Determine number of SNPs from length of risk_alleles vector
  n_snps <- length(risk_alleles)
  
  # Empty list to store results for each SNP
  tdt_results <- vector("list", n_snps)
  
  # Loop through each SNP and run compute_tdt() for that SNP + store the result
  for (i in 1:n_snps) {
    tdt_results[[i]] <- compute_tdt(trios, snp = i, risk_allele = risk_alleles[i])
  }
  
  # Label results as 'SNP1', 'SNP2' etc
  names(tdt_results) <- paste0("SNP", 1:n_snps)
  
  # Return results
  tdt_results
}


#####################################################
# 7. Function that performs the simulations and computes the TDT Z-scores

simulate_tdt_ld_correlation <- function(
  n_sim = 1000, # number of replicate trios
  n_reps = 2, # number of Monte Carlo simulations
  maf1 = 0.3, # maf1, maf2: minor allele frequencies for loci 1 and 2
  maf2 = 0.3, 
  r2 = 0.5, # linkage disequilibrium between the two loci
  p_transmit = 0.5 # probability of transmitting 'A' allele at SNP1 (0.5 = Mendelian)
) {
  # Empty matrix to store results
  z_scores <- matrix(NA, nrow = n_reps, ncol = 2)

  for (rep in 1:n_reps) {
    # Simulate trios
    trios <- simulate_trios(n_sim, maf1, maf2, r2, p_transmit)
    # Compute TDT Z-scores
    tdt <- compute_tdt_multi(trios, risk_alleles = c("A","B"))
    z_scores[rep, 1] <- tdt$SNP1$Z
    z_scores[rep, 2] <- tdt$SNP2$Z
  }

  # Convert to data.frame
  z_scores_df <- as.data.frame(z_scores)
  colnames(z_scores_df) <- c("Z_SNP1", "Z_SNP2")

  # Compute correlation between Z-scores
  cor_val <- cor(z_scores_df$Z_SNP1, z_scores_df$Z_SNP2)

  list(
    z_scores = z_scores_df,
    correlation = cor_val,
    r2 = r2
  )
}

#####################################################
# 8. Perform the simulations

set.seed(123)

start_time <- proc.time()

# Run across various MAF, r2 and TRD

# Create a grid of combinations to test
param_grid <- expand.grid(
  maf1 = c(0.01, 0.05, 0.1, 0.3, 0.5),
  maf2 = c(0.01, 0.05, 0.1, 0.3, 0.5),
  r2 = seq(0, 1, by = 0.1),
  p_transmit = c(0.5, 0.51, 0.6, 0.7, 0.8)
)

# Store results
results_df <- data.frame(
  maf1 = numeric(),
  maf2 = numeric(),
  r2 = numeric(),
  p_transmit = numeric(),
  correlation = numeric()
)

# Run simulations for each combination of parameters
for (i in 1:nrow(param_grid)) {
  # Try to run simulation, skip if r2 not achievable
  # Sometimes the MAF combinations make certain r2 values mathematically impossible
  res <- tryCatch({
    simulate_tdt_ld_correlation(
      n_sim = 1000,
      n_reps = 1000,
      maf1 = param_grid$maf1[i],
      maf2 = param_grid$maf2[i],
      r2 = param_grid$r2[i],
      p_transmit = param_grid$p_transmit[i]
    )
  }, error = function(e) {
    return(NULL)
  })
  
  # Only add result if simulation succeeded
  if (!is.null(res)) {
    results_df <- rbind(results_df, data.frame(
      maf1 = param_grid$maf1[i],
      maf2 = param_grid$maf2[i],
      r2 = param_grid$r2[i],
      p_transmit = param_grid$p_transmit[i],
      correlation = res$correlation
    ))
  } else {
    cat("Skipping maf1=", param_grid$maf1[i], 
        ", maf2=", param_grid$maf2[i], 
        ", r2=", param_grid$r2[i], 
        ", p_transmit=", param_grid$p_transmit[i],
        " (not achievable)\n", sep="")
  }
}

#####################################################
# 9. Create plots

# Compare the correlations (r against TDT Z-score correlation)
results_df$r <- sqrt(results_df$r2)

# Add omega column for clearer labeling
# ω (omega) quantifies deviation from expected 0.5 segregation ratio (p_transmit = 0.5 + ω)
results_df$omega <- results_df$p_transmit - 0.5

# Save results
write.table(results_df, "simulation_results.txt",
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Plot r vs TDT correlation faceted by both MAF1 and omega
p1 <- ggplot(results_df, aes(x = r, y = correlation, color = factor(maf2))) +
  geom_point(size = 2) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_grid(omega ~ maf1, 
             labeller = labeller(
               maf1 = as_labeller(function(x) paste0("MAF[SNP1]~`=`~", x), label_parsed),
               omega = as_labeller(function(x) paste0("omega~`=`~", x), label_parsed)
             )) +
  # Remove trailing zeros on x and y labels
  scale_x_continuous(labels = function(x) format(x, nsmall = 0, drop0trailing = TRUE)) +
  scale_y_continuous(labels = function(x) format(x, nsmall = 0, drop0trailing = TRUE)) +
  # Custom colour palette
  scale_color_manual(values = c("#a6d854", "#e78ac3", "#8da0cb", "#fc8d62", "#66c2a5")) +
  labs(
    x = "LD between SNP pairs (r)",
    y = "TDT Z-Score Correlation",
    color = expression(MAF[SNP2])
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

# Save plot
ggsave("tdt_ld_correlation.png", p1, width = 10, height = 10)

#####################################################
# 10. Overall correlation between r and TDT correlation
overall_cor <- cor(results_df$r, results_df$correlation)
cat("\nOverall correlation between r and TDT correlation:", round(overall_cor, 3), "\n")

end_time <- proc.time()
cat("\nTotal runtime:", round((end_time - start_time)[3] / 60, 2), "minutes\n")
