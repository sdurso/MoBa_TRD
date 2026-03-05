#!/usr/bin/env Rscript

# R script that computes and plots statistical power curves 
# for the Transmission Disequilibrium Test (TDT)
# under varying heterozygote frequencies and transmission ratio distortion (TRD; omega)

# Output:
#   e.g. Power_curve_moba_n_43841.png

###########################################################################################

library("reshape2")
library("ggplot2")

###########################################################################################
# Heterozygote frequencies
het_freq_list <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Transmission distortion parameter (omega)
w_list <- seq(0.01, 0.05, by = 0.0005)

# Set N a the number of trios - 1
n = 43840 

# Computes the critical chi-square value corresponding to the genome-wide significance threshold
null = qchisq(p= 0.00000005, df=1, ncp= 0, lower.tail = F, log.p=F)

# Set up results table
tdt_power <- matrix(data = NA, nrow = length(w_list), ncol=5)
cnames <- c("10%", "20%", "30%", "40%", "50%")
colnames(tdt_power) = cnames
rownames(tdt_power) = w_list

for (k in seq(length(het_freq_list))) {
  for (l in seq(length(w_list))) {
    het_freq <- het_freq_list[k]
    w <- w_list[l]
    # Non-centrality parameter for TDT
    ncp_alt = 4*(n*2*het_freq)*w^2
    power = pchisq(q = null , df=1, ncp=ncp_alt, lower.tail= F, log.p=F)
    tdt_power[l,k] <- power
  }
}

# Reshape for plotting
mdf <- melt(tdt_power, value.name="Power", varnames=c("TDT", "Het_Freq"))

# Plot
png(file="Power_curve_moba_n_43841.png", width=900, height=700,res=200)
ggplot(data=mdf, aes(x=TDT, y=Power, group = Het_Freq, colour = Het_Freq)) +
  labs(y ="Statistical Power", 
       x=expression(paste("Amount of Transmission Distortion (", omega,")"))) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "grey") +
  geom_line(size=1) +
  theme_bw() + 
  scale_colour_brewer(palette = "Set2") + # colour blind friendly
  guides(col = guide_legend(title = "Heterozygote \n Frequency")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        text=element_text(size=12),
        legend.title=element_text(size=10),
        panel.grid.minor = element_blank(), 
        axis.title.y = element_text(vjust = +2),
        axis.title.x = element_text(vjust = -0.75),
        axis.line = element_line(colour = "black"))
dev.off()
