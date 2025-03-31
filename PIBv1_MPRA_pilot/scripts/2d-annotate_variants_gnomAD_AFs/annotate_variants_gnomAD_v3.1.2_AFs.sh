#!/bin/sh
#SBATCH -a 1-22
#SBATCH --time=6:00:00
##SBATCH --partition=transfer
#SBATCH --mem=64G --cpus-per-task=1

# set i to job index
j=${SLURM_ARRAY_TASK_ID}

# run main script
# bgzip -dk ../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed.gz
tabix -h /home/sr2446/scratch60/gnomAD.genomes.v3.1.2.temp/gnomad.genomes.v3.1.2.sites.chr${j}.vcf.bgz -R ../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed >| ../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_lift38_chr-gnomad.genomes.v3.1.2.sites.chr${j}.vcf
bgzip ../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_lift38_chr-gnomad.genomes.v3.1.2.sites.chr${j}.vcf
# rm ../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed
