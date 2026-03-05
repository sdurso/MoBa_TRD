#!/bin/bash
#SBATCH --job-name=plink_tdt_example
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

# -----------------------------
# Example PLINK TDT pipeline
# -----------------------------
# This script shows how to run TDT and POO-TDT using PLINK
# Usage:
# sbatch plink_tdt.sh
# Adjust paths, file names, and SLURM settings for your own system
# -----------------------------

module load plink

# Set working directory (example placeholder)
WORKDIR="/path/to/your/data"
cd "$WORKDIR" || exit 1

# -----------------------------
# Example input files
# -----------------------------
BFILE="bfile"                       # PLINK binary fileset (bfile.bed/bim/fam)
KEEP_ALL="keep_trios.txt"           # List of trios to keep
KEEP_FEMALE="keep_female_trios.txt" # List of trios with female offspring
KEEP_MALE="keep_male_trios.txt"     # List of trios with male offspring
PHENO_FILE="pheno.txt"              # Phenotype file (to be generated)

# -----------------------------
# Create phenotype file
# -----------------------------
# Code all individuals as affected (2)
# Format: FID IID PHENOTYPE
awk '{print $1, $2, 2}' "${BFILE}.fam" > "$PHENO_FILE"

# -----------------------------
# Run TDT analyses
# -----------------------------

# Combined TDT (all offspring)
plink --bfile "$BFILE" \
      --pheno "$PHENO_FILE" \
      --keep "$KEEP_ALL" \
      --tdt \
      --out tdt_all_off

# Offspring sex–stratified TDT
plink --bfile "$BFILE" \
      --pheno "$PHENO_FILE" \
      --keep "$KEEP_FEMALE" \
      --tdt \
      --out tdt_female_off

plink --bfile "$BFILE" \
      --pheno "$PHENO_FILE" \
      --keep "$KEEP_MALE" \
      --tdt \
      --out tdt_male_off

# Parent-of-origin TDT
plink --bfile "$BFILE" \
      --pheno "$PHENO_FILE" \
      --keep "$KEEP_ALL" \
      --tdt \
      --poo \
      --out tdt_parent_of_origin

