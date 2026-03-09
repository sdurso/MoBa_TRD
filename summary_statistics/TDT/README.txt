================================================================================
TDT SUMMARY STATISTICS
================================================================================

Summary statistics from four TDT analyses in MoBa:
  - Combined (all offspring, autosomes + X)
  - Male-stratified (male offspring only, autosomes + X)
  - Female-stratified (female offspring only, autosomes + X)
  - Parent-of-origin (all offspring, autosomes only)

================================================================================
FILES
================================================================================

Summary statistics prior to DENTIST filtering:
  tdt_all_off_sumstats_DENTIST.txt
  tdt_male_off_sumstats_DENTIST.txt
  tdt_female_off_sumstats_DENTIST.txt
  tdt_parent_of_origin_sumstats_DENTIST.txt

Summary statistics after DENTIST filtering (retaining DENTIST_P > 5e-8):
  tdt_all_off_sumstats_DENTIST_pass.txt
  tdt_male_off_sumstats_DENTIST_pass.txt
  tdt_female_off_sumstats_DENTIST_pass.txt
  tdt_parent_of_origin_sumstats_DENTIST_pass.txt

Prefixes:
  tdt_all_off_           Combined analysis
  tdt_male_off_          Male-stratified analysis
  tdt_female_off_        Female-stratified analysis
  tdt_parent_of_origin_  Parent-of-origin analysis

================================================================================
COLUMN DESCRIPTIONS
================================================================================

Combined, Male-Stratified, and Female-Stratified Analyses:
SNP         SNP identifier
CHR         Chromosome
BP          Base pair position (GRCh37)
A1          Minor allele in MoBa
A2          Major allele in MoBa
MAF         Minor allele frequency in MoBa
T           Transmission count of the minor allele
U           Transmission count of the major allele
Total       Total tranmsission count (T + U)
OR          TDT odds (T/U)
CHISQ       TDT chi-squared test statistic
P           TDT p-value
DENTIST_P   DENTIST p-value
HWE_P       Hardy-Weinberg equilibrium p-value (founder genotypes)
FDR         Benjamini-Hochberg FDR adjusted p-value, calculated post-DENTIST filtering 

Parent-of-Origin Analysis:
SNP             SNP identifier
CHR             Chromosome
BP              Base pair position (GRCh37)
A1              Minor allele in MoBa
A2              Major allele in MoBa
MAF             Minor allele frequency in MoBa
minor_MAT_count Maternal transmission count of minor allele
major_MAT_count Maternal transmission count of major allele
CHISQ_MAT       Maternal TDT chi-squared statistic
MAT_OR          Maternal odds (minor_MAT_count/major_MAT_count)
P_MAT           Maternal TDT p-value
minor_PAT_count Paternal transmission count of minor allele
major_PAT_count Paternal transmission count of major allele
CHISQ_PAT       Paternal TDT chi-squared statistic
PAT_OR          Paternal odds (minor_PAT_count/major_PAT_count)
P_PAT           Paternal TDT p-value
Total           Total transmission count (minor_MAT_count + major_MAT_count + minor_PAT_count + major_PAT_count)
Z_POO           Parent-of-origin z-score (maternal vs paternal)
P_POO           Parent-of-origin p-value
DENTIST_P       DENTIST p-value
HWE_P           Hardy-Weinberg equilibrium p-value (founder genotypes)
FDR             Benjamini-Hochberg FDR adjusted p-value, calculated post-DENTIST filtering 

================================================================================
