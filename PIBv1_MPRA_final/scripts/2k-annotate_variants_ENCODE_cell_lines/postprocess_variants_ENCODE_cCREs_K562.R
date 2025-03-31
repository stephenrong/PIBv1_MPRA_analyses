#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)

# load VCF 
introgressed_variants_ENCODE_cCREs_K562 <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)

# load BED
ENCODE_cCREs_K562_tb <- as_tibble(fread("../../../../Datasets/gene_regulation_element_catalogs/ENCODE_SCREEN_V4_cCREs/data_cleanup/sample_specific/ENCODE_SCREEN_V4_cCREs_cell_lines_lift37/K562_ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.lift37.nochr.sort.bed.gz"))
ENCODE_cCREs_K562_tb <- ENCODE_cCREs_K562_tb[,c(1,2,3,10)]
ENCODE_cCREs_K562_tb <- unique(ENCODE_cCREs_K562_tb)
names(ENCODE_cCREs_K562_tb) <- c("seqnames", "start", "end", "ENCODE_cCREs_K562_orig_class")
ENCODE_cCREs_K562_tb <- ENCODE_cCREs_K562_tb %>% 
	mutate(start = start + 1) %>%  # translate 0-indexed to 1-indexed
	mutate(ENCODE_cCREs_K562_orig_element = paste(paste(seqnames, start, end, sep="_"), "K562", ENCODE_cCREs_K562_orig_class, sep="-"))

# overlap
overlap_temp <- as_tibble(findOverlaps(GRanges(introgressed_variants_ENCODE_cCREs_K562), GRanges(ENCODE_cCREs_K562_tb)))
overlap_temp$VariantID <- introgressed_variants_ENCODE_cCREs_K562$VariantID[overlap_temp$queryHits]
overlap_temp$ENCODE_cCREs_K562_orig_class <- ENCODE_cCREs_K562_tb$ENCODE_cCREs_K562_orig_class[overlap_temp$subjectHits]
overlap_temp$ENCODE_cCREs_K562_orig_element <- ENCODE_cCREs_K562_tb$ENCODE_cCREs_K562_orig_element[overlap_temp$subjectHits]
overlap_temp <- overlap_temp %>% 
	dplyr::select(-queryHits, -subjectHits)
introgressed_variants_ENCODE_cCREs_K562 <- introgressed_variants_ENCODE_cCREs_K562 %>% 
	left_join(overlap_temp)

# save names
introgressed_variants_ENCODE_cCREs_K562_names <- names(introgressed_variants_ENCODE_cCREs_K562)[11:length(introgressed_variants_ENCODE_cCREs_K562)]

# collapse multiple entries
introgressed_variants_ENCODE_cCREs_K562 <- introgressed_variants_ENCODE_cCREs_K562 %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(introgressed_variants_ENCODE_cCREs_K562_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup()

# add summary cols
introgressed_variants_ENCODE_cCREs_K562 <- introgressed_variants_ENCODE_cCREs_K562 %>% 
	mutate(ENCODE_cCREs_K562_summ_Category = ENCODE_cCREs_K562_orig_class) %>% 
	mutate(ENCODE_cCREs_K562_summ_Description = ifelse(ENCODE_cCREs_K562_orig_class %in% c("Low-DNase"), "NA", ENCODE_cCREs_K562_orig_element)) %>% 
	mutate(ENCODE_cCREs_K562_summ_Summary = ifelse(ENCODE_cCREs_K562_orig_class %in% c("Low-DNase"), "NA", ENCODE_cCREs_K562_orig_class))

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_inter <- intersect(names(introgressed_variants_tb), names(introgressed_variants_ENCODE_cCREs_K562))
introgressed_variants_ENCODE_cCREs_K562 <- introgressed_variants_tb %>% mutate_at(names_inter, as.character) %>% 
	left_join(introgressed_variants_ENCODE_cCREs_K562 %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(introgressed_variants_ENCODE_cCREs_K562, gzfile("../../results/2k-annotate_variants_ENCODE_cell_lines/introgressed_variants_ENCODE_cCREs_K562.txt.gz"))

# save summ variants
introgressed_variants_ENCODE_cCREs_K562_summ <- introgressed_variants_ENCODE_cCREs_K562 %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_ENCODE_cCREs_K562_summ, gzfile("../../results/2k-annotate_variants_ENCODE_cell_lines/introgressed_variants_ENCODE_cCREs_K562_summ.txt.gz"))
sort(table(introgressed_variants_ENCODE_cCREs_K562_summ$ENCODE_cCREs_K562_summ_Summary))
