#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load VCF
introgressed_variants_ClinVar <- read.vcfR("../../results/2h-annotate_variants_ClinVar_categories/introgressed_variants-clinvar_20230326.vcf.gz")

# split INFO
introgressed_variants_ClinVar <- vcfR2tidy(introgressed_variants_ClinVar, info_only=TRUE)$fix

# filter SNPs
introgressed_variants_ClinVar <- introgressed_variants_ClinVar %>% 
	filter(REF %in% c("A", "C", "G", "T"), ALT %in% c("A", "C", "G", "T"))

# save names
names(introgressed_variants_ClinVar)[c(3,6:length(introgressed_variants_ClinVar))] <- paste("ClinVar_20230326_orig", names(introgressed_variants_ClinVar)[c(3,6:length(introgressed_variants_ClinVar))], sep="_")
introgressed_variants_ClinVar_names <- names(introgressed_variants_ClinVar)[c(3,6:length(introgressed_variants_ClinVar))]

# add GRange and Variant cols
introgressed_variants_ClinVar <- introgressed_variants_ClinVar %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, all_of(introgressed_variants_ClinVar_names))

# add summary cols
introgressed_variants_ClinVar <- introgressed_variants_ClinVar %>% 
	mutate(ClinVar_20230326_summ_Summary = ifelse(grepl("^criteria_provided", ClinVar_20230326_orig_CLNREVSTAT), ifelse(grepl("Pathogenic|Likely_pathogenic", ClinVar_20230326_orig_CLNSIG), "Pathogenic/LP", ifelse(grepl("Benign|Likely_benign", ClinVar_20230326_orig_CLNSIG), "Benign/LB", "VUS/NCP/Other")), "VUS/NCP/Other"))

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_ClinVar))
introgressed_variants_ClinVar <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_ClinVar %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(introgressed_variants_ClinVar, gzfile("../../results/2h-annotate_variants_ClinVar_categories/introgressed_variants_ClinVar.txt.gz"))

# save summ variants
introgressed_variants_ClinVar_summ <- introgressed_variants_ClinVar %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_ClinVar_summ, gzfile("../../results/2h-annotate_variants_ClinVar_categories/introgressed_variants_ClinVar_summ.txt.gz"))
sort(table(introgressed_variants_ClinVar_summ$ClinVar_20230326_summ_Summary))
