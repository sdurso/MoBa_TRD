#!/bin/bash
#SBATCH --job-name=dentist_example
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00

# -----------------------------
# Example DENTIST pipeline
# -----------------------------
# This script shows how to run DENTIST QC of TDT summary statistics

# Usage:
#for chr in {1..23}; do
#    for tdt in dentist_tdt_all_off_chr${chr}.txt \
#               dentist_tdt_male_off_chr${chr}.txt \
#               dentist_tdt_female_off_chr${chr}.txt; do
#        sbatch dentist.sh ${tdt} ${chr}
#    done
#done

#for chr in {1..22}; do
#    for tdt in dentist_tdt_parent_of_origin_chr${chr}.txt; do
#        sbatch dentist.sh ${tdt} ${chr}
#    done
#done

# Adjust paths, file names, and SLURM settings for your own system

# Output files *.full.txt
1_dentist_tdt_all_off_chr1.txt.DENTIST.full.txt


# -----------------------------
# Example input files
# -----------------------------
BFILE="bfile" # PLINK binary fileset (bfile.bed/bim/fam)

module load DENTIST

# Set working directory (example placeholder)
WORKDIR="/path/to/your/data"
cd "$WORKDIR" || exit 1

# Arguments passed to script
TDT_FILE=${1}   # e.g. dentist_tdt_all_off_chr1.txt
CHR=${2}        # e.g. 1

# Run DENTIST
DENTIST --gwas-summary ${TDT_FILE} \
        --bfile ${BFILE} \
        --chrID ${CHR} \
        --thread-num 20 \
        --out ${CHR}_${TDT_FILE%.txt}
