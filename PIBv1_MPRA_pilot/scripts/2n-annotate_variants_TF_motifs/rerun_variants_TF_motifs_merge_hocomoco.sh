#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=32G --cpus-per-task=1

module load R/4.2.0-foss-2020b
Rscript --vanilla rerun_variants_TF_motifs_merge_hocomoco.R
