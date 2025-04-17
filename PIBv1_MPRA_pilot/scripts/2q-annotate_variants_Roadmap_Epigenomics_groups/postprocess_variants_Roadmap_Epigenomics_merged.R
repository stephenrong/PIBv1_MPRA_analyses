#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(rtracklayer)
source("../shared_functions/seqinfo_fix_change.R")

# load variants
adaptive_variants_tb <- as_tibble(readRDS("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.rds"))[,1:10]

# load Roadmap_Epigenomics_merged
Roadmap_Epigenomics_merged_tb <- as_tibble(fread(paste0("../../../Datasets/gene_regulation_element_catalogs/Roadmap_Epigenomics_elements/data_cleanup/GRCh37_merged/Roadmap_Epigenomics_merged2_groups.txt.gz")))

# rename group
Roadmap_Epigenomics_merged_tb <- Roadmap_Epigenomics_merged_tb %>% 
	mutate(group = gsub(" ", "_", group)) %>% 
	mutate(group = gsub("\\.", "", group)) %>% 
	mutate(group = gsub("&", "and", group))

# Roadmap_Epigenomics_merged names
names(Roadmap_Epigenomics_merged_tb)[6:length(Roadmap_Epigenomics_merged_tb)] <- paste("Roadmap_Epigenomics_merged_orig", names(Roadmap_Epigenomics_merged_tb)[6:length(Roadmap_Epigenomics_merged_tb)], sep="_")
adaptive_variants_Roadmap_Epigenomics_merged_names <- names(Roadmap_Epigenomics_merged_tb)[6:length(Roadmap_Epigenomics_merged_tb)]

# overlap variants and Roadmap_Epigenomics_merged
adaptive_variants_Roadmap_Epigenomics_merged <- as_tibble(find_overlaps(GRanges(adaptive_variants_tb), GRanges(Roadmap_Epigenomics_merged_tb)))

# Roadmap_Epigenomics_merged names
adaptive_variants_Roadmap_Epigenomics_merged_names <- names(adaptive_variants_Roadmap_Epigenomics_merged)[11:length(adaptive_variants_Roadmap_Epigenomics_merged)]

# collapse by group
adaptive_variants_Roadmap_Epigenomics_merged_collapse <- adaptive_variants_Roadmap_Epigenomics_merged %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	dplyr::summarise_at(adaptive_variants_Roadmap_Epigenomics_merged_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup() %>% 
	mutate(Roadmap_Epigenomics_merged_summ_Group = Roadmap_Epigenomics_merged_orig_group) %>% 
	mutate(Roadmap_Epigenomics_merged_summ_Summary = "Roadmap_Epigenomics_merged_overlap") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, 
		starts_with("Roadmap_Epigenomics_merged_orig"), Roadmap_Epigenomics_merged_summ_Summary, Roadmap_Epigenomics_merged_summ_Group)

# pivot by group
adaptive_variants_Roadmap_Epigenomics_merged_pivot <- adaptive_variants_Roadmap_Epigenomics_merged %>% 
	mutate(Roadmap_Epigenomics_merged_orig_group_temp = Roadmap_Epigenomics_merged_orig_group) %>% 
	pivot_wider(names_from=Roadmap_Epigenomics_merged_orig_group, names_sep="-", names_prefix="Roadmap_Epigenomics_merged_summ_Group-", values_from=Roadmap_Epigenomics_merged_orig_group_temp)

# join together
adaptive_variants_Roadmap_Epigenomics_merged <- full_join(adaptive_variants_Roadmap_Epigenomics_merged_collapse, adaptive_variants_Roadmap_Epigenomics_merged_pivot)

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_Roadmap_Epigenomics_merged))
adaptive_variants_Roadmap_Epigenomics_merged <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_Roadmap_Epigenomics_merged %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_Roadmap_Epigenomics_merged, gzfile("../../results/2q-annotate_variants_Roadmap_Epigenomics_groups/adaptive_variants_Roadmap_Epigenomics_merged.txt.gz"))

# save summ variants
adaptive_variants_Roadmap_Epigenomics_merged_summ <- adaptive_variants_Roadmap_Epigenomics_merged %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_Roadmap_Epigenomics_merged_summ, gzfile("../../results/2q-annotate_variants_Roadmap_Epigenomics_groups/adaptive_variants_Roadmap_Epigenomics_merged_summ.txt.gz"))
sort(table(adaptive_variants_Roadmap_Epigenomics_merged_summ$Roadmap_Epigenomics_merged_summ_Summary))
