#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --partition=week
#SBATCH --mem=256G --cpus-per-task=1

module load R/4.2.0-foss-2020b
Rscript --vanilla rerun_variants_TF_motifs_merge_jaspar.R
