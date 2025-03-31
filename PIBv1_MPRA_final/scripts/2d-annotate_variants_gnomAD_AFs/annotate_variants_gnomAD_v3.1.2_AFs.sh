#!/bin/bash
#SBATCH -a 1-22
#SBATCH --time=24:00:00
##SBATCH --partition=transfer
#SBATCH --mem=64G --cpus-per-task=1

# set i to job index
j=${SLURM_ARRAY_TASK_ID}

# run main script
tabix -h /vast/palmer/scratch/reilly/sr2446/gnomad_temp/gnomad.genomes.v3.1.2.sites.chr${j}.vcf.bgz -R ../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/introgressed_variants_lift38_chr.bed >| ../../results/2d-annotate_variants_gnomAD_AFs/introgressed_variants_lift38_chr-gnomad.genomes.v3.1.2.sites.chr${j}.vcf
bgzip -f ../../results/2d-annotate_variants_gnomAD_AFs/introgressed_variants_lift38_chr-gnomad.genomes.v3.1.2.sites.chr${j}.vcf
