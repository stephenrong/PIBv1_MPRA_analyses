#!/bin/sh
#SBATCH -a 1-4
#SBATCH --time=12:00:00
#SBATCH --mem=32G --cpus-per-task=1

# set i to job index
j=${SLURM_ARRAY_TASK_ID}

# run main script
if [ "$j" == "1" ]
then
	tabix -h ../../../../Datasets/human_popgen_genomes/archaic_high_cov_genomes/altai_denisovan/altai_denisovan_masked_norm.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.bed >| ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-altai_denisovan_masked_norm.vcf
	bgzip -f ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-altai_denisovan_masked_norm.vcf
fi

if [ "$j" == "2" ]
then
	tabix -h ../../../../Datasets/human_popgen_genomes/archaic_high_cov_genomes/altai_neanderthal/altai_neanderthal_masked_norm.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.bed >| ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-altai_neanderthal_masked_norm.vcf
	bgzip -f ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-altai_neanderthal_masked_norm.vcf
fi

if [ "$j" == "3" ]
then
	tabix -h ../../../../Datasets/human_popgen_genomes/archaic_high_cov_genomes/vindija_neanderthal/vindija_neanderthal_masked_norm.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.bed >| ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-vindija_neanderthal_masked_norm.vcf
	bgzip -f ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-vindija_neanderthal_masked_norm.vcf
fi

if [ "$j" == "4" ]
then
	tabix -h ../../../../Datasets/human_popgen_genomes/archaic_high_cov_genomes/chagyrskaya_neanderthal/chagyrskaya_neanderthal_masked_norm.vcf.gz -R ../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.bed >| ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-chagyrskaya_neanderthal_masked_norm.vcf
	bgzip -f ../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-chagyrskaya_neanderthal_masked_norm.vcf
fi
