#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load VCF
adaptive_variants_1KGP_phase3_AFs <- read.vcfR("../../results/2c-annotate_variants_1KGP_AFs/adaptive_variants-ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz")

# split INFO
adaptive_variants_1KGP_phase3_AFs <- vcfR2tidy(adaptive_variants_1KGP_phase3_AFs, info_only=TRUE)$fix

# filter SNPs
adaptive_variants_1KGP_phase3_AFs <- adaptive_variants_1KGP_phase3_AFs %>% 
	filter(REF %in% c("A", "C", "G", "T"), ALT %in% c("A", "C", "G", "T"))

# add prefix and save names
names(adaptive_variants_1KGP_phase3_AFs)[c(3,6:length(adaptive_variants_1KGP_phase3_AFs))] <- paste("1KGP_phase3_orig", names(adaptive_variants_1KGP_phase3_AFs)[c(3,6:length(adaptive_variants_1KGP_phase3_AFs))], sep="_")
adaptive_variants_1KGP_phase3_AFs_names <- names(adaptive_variants_1KGP_phase3_AFs)[c(3,6:length(adaptive_variants_1KGP_phase3_AFs))]

# add GRange and Variant cols
adaptive_variants_1KGP_phase3_AFs <- adaptive_variants_1KGP_phase3_AFs %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, all_of(adaptive_variants_1KGP_phase3_AFs_names))

# add summary AF, AC, AN cols
adaptive_variants_1KGP_phase3_AFs <- adaptive_variants_1KGP_phase3_AFs %>% 
	mutate(
		`1KGP_phase3_summ_AF` = as.numeric(`1KGP_phase3_orig_AF`), 
		`1KGP_phase3_summ_AC` = as.numeric(`1KGP_phase3_orig_AC`), 
		`1KGP_phase3_summ_AN` = as.numeric(`1KGP_phase3_orig_AN`), 
		`1KGP_phase3_summ_EAS_AF` = as.numeric(`1KGP_phase3_orig_EAS_AF`), 
		`1KGP_phase3_summ_EAS_AC` = round(as.numeric(`1KGP_phase3_orig_EAS_AF`)*2*504), 
		`1KGP_phase3_summ_EAS_AN` = 2*504, 
		`1KGP_phase3_summ_EUR_AF` = as.numeric(`1KGP_phase3_orig_EUR_AF`), 
		`1KGP_phase3_summ_EUR_AC` = round(as.numeric(`1KGP_phase3_orig_EUR_AF`)*2*503), 
		`1KGP_phase3_summ_EUR_AN` = 2*503, 
		`1KGP_phase3_summ_AFR_AF` = as.numeric(`1KGP_phase3_orig_AFR_AF`), 
		`1KGP_phase3_summ_AFR_AC` = round(as.numeric(`1KGP_phase3_orig_AFR_AF`)*2*661), 
		`1KGP_phase3_summ_AFR_AN` = 2*661, 
		`1KGP_phase3_summ_AMR_AF` = as.numeric(`1KGP_phase3_orig_AMR_AF`), 
		`1KGP_phase3_summ_AMR_AC` = round(as.numeric(`1KGP_phase3_orig_AMR_AF`)*2*347), 
		`1KGP_phase3_summ_AMR_AN` = 2*347, 
		`1KGP_phase3_summ_SAS_AF` = as.numeric(`1KGP_phase3_orig_SAS_AF`), 
		`1KGP_phase3_summ_SAS_AC` = round(as.numeric(`1KGP_phase3_orig_SAS_AF`)*2*489), 
		`1KGP_phase3_summ_SAS_AN` = 2*489
	) %>% 
	mutate(
		`1KGP_phase3_summ_AF` = round(`1KGP_phase3_summ_AC`/`1KGP_phase3_summ_AN`, 6), 
		`1KGP_phase3_summ_EAS_AF` = round(`1KGP_phase3_summ_EAS_AC`/`1KGP_phase3_summ_EAS_AN`, 6), 
		`1KGP_phase3_summ_EUR_AF` = round(`1KGP_phase3_summ_EUR_AC`/`1KGP_phase3_summ_EUR_AN`, 6), 
		`1KGP_phase3_summ_AFR_AF` = round(`1KGP_phase3_summ_AFR_AC`/`1KGP_phase3_summ_AFR_AN`, 6), 
		`1KGP_phase3_summ_AMR_AF` = round(`1KGP_phase3_summ_AMR_AC`/`1KGP_phase3_summ_AMR_AN`, 6), 
		`1KGP_phase3_summ_SAS_AF` = round(`1KGP_phase3_summ_SAS_AC`/`1KGP_phase3_summ_SAS_AN`, 6)
	)

# save only summary cols
adaptive_variants_1KGP_phase3_AFs <- adaptive_variants_1KGP_phase3_AFs %>% 
	dplyr::select(-starts_with("1KGP_phase3_orig_"))

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_1KGP_phase3_AFs))
adaptive_variants_1KGP_phase3_AFs <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_1KGP_phase3_AFs %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_1KGP_phase3_AFs, gzfile("../../results/2c-annotate_variants_1KGP_AFs/adaptive_variants_1KGP_phase3_AFs.txt.gz"))

# save summ variants
adaptive_variants_1KGP_phase3_AFs_summ <- adaptive_variants_1KGP_phase3_AFs %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_1KGP_phase3_AFs_summ, gzfile("../../results/2c-annotate_variants_1KGP_AFs/adaptive_variants_1KGP_phase3_AFs_summ.txt.gz"))
