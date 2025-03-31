#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load TXT
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)

UKBB_94traits_finemap_pip0.1 <- as_tibble(fread("../../../../Datasets/fine_mapped_regions_gwas/Kanai_et_2021_fine_mapped_ukbb/data_cleanup/UKBB_94traits_finemap_pip0.1.txt.gz"))
names(UKBB_94traits_finemap_pip0.1) <- gsub("UKBB", "UKBB_finemap_orig", names(UKBB_94traits_finemap_pip0.1))
write_tsv(UKBB_94traits_finemap_pip0.1, gzfile("../../results/2m-annotate_variants_UKBB_BBJ_finemap/UKBB_94traits_finemap_pip0.1.txt.gz"))

# filter variants
adaptive_variants_UKBB_finemap <- UKBB_94traits_finemap_pip0.1 %>% 
	filter(VariantID %in% adaptive_variants_tb$VariantID)

# add temp col
adaptive_variants_UKBB_finemap <- adaptive_variants_UKBB_finemap %>% 
	mutate(UKBB_finemap_summ_Description = paste(UKBB_finemap_orig_trait, UKBB_finemap_orig_method, UKBB_finemap_orig_region, UKBB_finemap_orig_cs_id, paste0("PIP=", UKBB_finemap_orig_pip), sep=","))

# save names
adaptive_variants_UKBB_finemap_names <- names(adaptive_variants_UKBB_finemap)[11:length(adaptive_variants_UKBB_finemap)]

# collapse multiple entries
adaptive_variants_UKBB_finemap <- adaptive_variants_UKBB_finemap %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	arrange(desc(UKBB_finemap_orig_pip)) %>% 
	dplyr::summarise_at(adaptive_variants_UKBB_finemap_names, function(x) {paste(x, collapse=";")}) %>% 
	ungroup()

# add summary cols
adaptive_variants_UKBB_finemap <- adaptive_variants_UKBB_finemap %>% 
	rowwise() %>% 
	mutate(UKBB_finemap_summ_Summary = paste(unique(strsplit(UKBB_finemap_orig_trait, split=",")[[1]]), collapse=",")) %>% 
	ungroup()

# merge to full variants
names_inter <- intersect(names(adaptive_variants_tb), names(adaptive_variants_UKBB_finemap))
adaptive_variants_UKBB_finemap <- adaptive_variants_tb %>% mutate_at(names_inter, as.character)%>% 
	left_join(adaptive_variants_UKBB_finemap %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(adaptive_variants_UKBB_finemap, gzfile("../../results/2m-annotate_variants_UKBB_BBJ_finemap/adaptive_variants_UKBB_finemap.txt.gz"))

# save summ variants
adaptive_variants_UKBB_finemap_summ <- adaptive_variants_UKBB_finemap %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_UKBB_finemap_summ, gzfile("../../results/2m-annotate_variants_UKBB_BBJ_finemap/adaptive_variants_UKBB_finemap_summ.txt.gz"))
sort(table(adaptive_variants_UKBB_finemap_summ$UKBB_finemap_summ_Summary))
