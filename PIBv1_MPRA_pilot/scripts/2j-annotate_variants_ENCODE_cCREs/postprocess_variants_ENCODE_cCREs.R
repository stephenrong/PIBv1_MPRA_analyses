#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
# library(vcfR)

# load VCF 
adaptive_variants_ENCODE_cCREs <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)

# load BED
ENCODE_cCREs_tb <- as_tibble(fread("../../../../Datasets/gene_regulation_element_catalogs/ENCODE_SCREEN_V4_cCREs/data_cleanup/sample_agnostic/ENCODE_SCREEN_V4_cCREs_lift37/GRCh38-cCREs.V4.lift37.nochr.sort.bed.gz"))
ENCODE_cCREs_tb <- ENCODE_cCREs_tb[,c(1,2,3,6)]
ENCODE_cCREs_tb <- unique(ENCODE_cCREs_tb)
names(ENCODE_cCREs_tb) <- c("seqnames", "start", "end", "ENCODE_cCREs_orig_class")
ENCODE_cCREs_tb <- ENCODE_cCREs_tb %>% 
	mutate(start = start + 1) %>%  # translate 0-indexed to 1-indexed
	mutate(ENCODE_cCREs_orig_element = paste(paste(seqnames, start, end, sep="_"), ENCODE_cCREs_orig_class, sep="-"))
write_tsv(ENCODE_cCREs_tb, gzfile("../../results/2j-annotate_variants_ENCODE_cCREs/ENCODE_cCREs_tb.txt.gz"))

# overlap
overlap_temp <- as_tibble(findOverlaps(GRanges(adaptive_variants_ENCODE_cCREs), GRanges(ENCODE_cCREs_gr)))
overlap_temp$VariantID <- adaptive_variants_ENCODE_cCREs$VariantID[overlap_temp$queryHits]
overlap_temp$ENCODE_cCREs_orig_class <- ENCODE_cCREs_gr$ENCODE_cCREs_orig_class[overlap_temp$subjectHits]
overlap_temp$ENCODE_cCREs_orig_element <- ENCODE_cCREs_gr$ENCODE_cCREs_orig_element[overlap_temp$subjectHits]
overlap_temp <- overlap_temp %>% 
	dplyr::select(-queryHits, -subjectHits)
adaptive_variants_ENCODE_cCREs <- adaptive_variants_ENCODE_cCREs %>% 
	left_join(overlap_temp)

# save names
adaptive_variants_ENCODE_cCREs_names <- names(adaptive_variants_ENCODE_cCREs)[11:length(adaptive_variants_ENCODE_cCREs)]

# collapse multiple entries
adaptive_variants_ENCODE_cCREs <- adaptive_variants_ENCODE_cCREs %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(adaptive_variants_ENCODE_cCREs_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup()

# add summary cols
adaptive_variants_ENCODE_cCREs <- adaptive_variants_ENCODE_cCREs %>% 
	mutate(ENCODE_cCREs_summ_Description = ENCODE_cCREs_orig_element) %>% 
	mutate(ENCODE_cCREs_summ_Summary = ENCODE_cCREs_orig_class)

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_inter <- intersect(names(adaptive_variants_tb), names(adaptive_variants_ENCODE_cCREs))
adaptive_variants_ENCODE_cCREs <- adaptive_variants_tb %>% mutate_at(names_inter, as.character) %>% 
	left_join(adaptive_variants_ENCODE_cCREs %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(adaptive_variants_ENCODE_cCREs, gzfile("../../results/2j-annotate_variants_ENCODE_cCREs/adaptive_variants_ENCODE_cCREs.txt.gz"))

# save summ variants
adaptive_variants_ENCODE_cCREs_summ <- adaptive_variants_ENCODE_cCREs %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_ENCODE_cCREs_summ, gzfile("../../results/2j-annotate_variants_ENCODE_cCREs/adaptive_variants_ENCODE_cCREs_summ.txt.gz"))
sort(table(adaptive_variants_ENCODE_cCREs_summ$ENCODE_cCREs_summ_Summary))
