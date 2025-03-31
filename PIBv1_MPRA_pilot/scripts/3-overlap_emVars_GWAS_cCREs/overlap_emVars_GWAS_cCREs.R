#!/bin/R

# load packages
library(tidyverse)
library(data.table)

# load bio packages
library(plyranges)

# load files
ENCODE_cCREs <- GRanges(as_tibble(fread("../../results/2j-annotate_variants_ENCODE_cCREs/ENCODE_cCREs_tb.txt.gz")))

UKBB_finemap_pip0.1 <- GRanges(as_tibble(fread("../../results/2m-annotate_variants_UKBB_BBJ_finemap/UKBB_94traits_finemap_pip0.1.txt.gz")))
names(mcols(UKBB_finemap_pip0.1))[1:5] <- paste0("UKBB_finemap_", names(mcols(UKBB_finemap_pip0.1))[1:5])

BBJ_finemap_pip0.1 <- GRanges(as_tibble(fread("../../results/2m-annotate_variants_UKBB_BBJ_finemap/BBJ_79traits_finemap_pip0.1.txt.gz")))
names(mcols(BBJ_finemap_pip0.1))[1:5] <- paste0("BBJ_finemap_", names(mcols(BBJ_finemap_pip0.1))[1:5])

mpra_ccre_variants <- GRanges(as_tibble(fread("../../results/3-merge_all_variant_annotations/mpra_ccre_variants_all_annotations_summ.txt.gz")))
names(mcols(mpra_ccre_variants))[1:5] <- paste0("MPRA_", names(mcols(mpra_ccre_variants))[1:5])

# get emVar 
mpra_emVars_variants <- mpra_ccre_variants %>% 
	filter(!is.na(MPRA_K562_summ_Summary) | !is.na(MPRA_Jurkat_summ_Summary))

# get cCRes + pip0.1
ENCODE_cCREs_UKBB_finemap_pip0.1 <- find_overlaps(ENCODE_cCREs, UKBB_finemap_pip0.1)

ENCODE_cCREs_BBJ_finemap_pip0.1 <- find_overlaps(ENCODE_cCREs, BBJ_finemap_pip0.1)

# get emVars + cCREs + pip0.1
mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1 <- find_overlaps(mpra_emVars_variants, ENCODE_cCREs_UKBB_finemap_pip0.1) %>% 
	as_tibble() %>% 
	arrange(-UKBB_finemap_orig_pip)

mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1 <- find_overlaps(mpra_emVars_variants, ENCODE_cCREs_BBJ_finemap_pip0.1) %>% 
	as_tibble() %>% 
	arrange(-BBJ_finemap_orig_pip)

# get abbreviated tables
mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_abbrev <- mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1 %>% 
	rename_with(~ gsub("X1KGP", "1KGP", .x)) %>% 
	dplyr::select(seqnames, start, end, width, strand, archaic_genotype_summ_Altai_Denisovan, archaic_genotype_summ_Altai_Neanderthal, archaic_genotype_summ_Vindija_Neanderthal, archaic_genotype_summ_Chagyrskaya_Neanderthal, `1KGP_phase3_summ_AF`, `1KGP_phase3_summ_AC`, `1KGP_phase3_summ_AN`, `gnomAD_v3.1.2_summ_AF`, `gnomAD_v3.1.2_summ_AC`, `gnomAD_v3.1.2_summ_AN`, MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_Jurkat_summ_Summary, nearest_gene_summ_Symbol, ABC_summ_Symbol, Roadmap_Epigenomics_summ_Symbol, EpiMap_summ_Symbol, ENCODE_cCREs_orig_class, ENCODE_cCREs_orig_element, UKBB_finemap_VariantID, UKBB_finemap_orig_method, UKBB_finemap_orig_trait, UKBB_finemap_orig_pip)

mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_abbrev <- mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1 %>% 
	rename_with(~ gsub("X1KGP", "1KGP", .x)) %>% 
	dplyr::select(seqnames, start, end, width, strand, archaic_genotype_summ_Altai_Denisovan, archaic_genotype_summ_Altai_Neanderthal, archaic_genotype_summ_Vindija_Neanderthal, archaic_genotype_summ_Chagyrskaya_Neanderthal, `1KGP_phase3_summ_AF`, `1KGP_phase3_summ_AC`, `1KGP_phase3_summ_AN`, `gnomAD_v3.1.2_summ_AF`, `gnomAD_v3.1.2_summ_AC`, `gnomAD_v3.1.2_summ_AN`, MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_Jurkat_summ_Summary, nearest_gene_summ_Symbol, ABC_summ_Symbol, Roadmap_Epigenomics_summ_Symbol, EpiMap_summ_Symbol, ENCODE_cCREs_orig_class, ENCODE_cCREs_orig_element, BBJ_finemap_VariantID, BBJ_finemap_orig_method, BBJ_finemap_orig_trait, BBJ_finemap_orig_pip)

# get emVars + pip0.1
mpra_emVars_variants_UKBB_finemap_pip0.1 <- find_overlaps(mpra_emVars_variants, UKBB_finemap_pip0.1) %>% 
	as_tibble() %>% 
	arrange(-UKBB_finemap_orig_pip)

mpra_emVars_variants_BBJ_finemap_pip0.1 <- find_overlaps(mpra_emVars_variants, BBJ_finemap_pip0.1) %>% 
	as_tibble() %>% 
	arrange(-BBJ_finemap_orig_pip)

# get abbreviated tables
mpra_emVars_variants_UKBB_finemap_pip0.1_abbrev <- mpra_emVars_variants_UKBB_finemap_pip0.1 %>% 
	rename_with(~ gsub("X1KGP", "1KGP", .x)) %>% 
	dplyr::select(seqnames, start, end, width, strand, archaic_genotype_summ_Altai_Denisovan, archaic_genotype_summ_Altai_Neanderthal, archaic_genotype_summ_Vindija_Neanderthal, archaic_genotype_summ_Chagyrskaya_Neanderthal, `1KGP_phase3_summ_AF`, `1KGP_phase3_summ_AC`, `1KGP_phase3_summ_AN`, `gnomAD_v3.1.2_summ_AF`, `gnomAD_v3.1.2_summ_AC`, `gnomAD_v3.1.2_summ_AN`, MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_Jurkat_summ_Summary, nearest_gene_summ_Symbol, ABC_summ_Symbol, Roadmap_Epigenomics_summ_Symbol, EpiMap_summ_Symbol, UKBB_finemap_VariantID, UKBB_finemap_orig_method, UKBB_finemap_orig_trait, UKBB_finemap_orig_pip)

mpra_emVars_variants_BBJ_finemap_pip0.1_abbrev <- mpra_emVars_variants_BBJ_finemap_pip0.1 %>% 
	rename_with(~ gsub("X1KGP", "1KGP", .x)) %>% 
	dplyr::select(seqnames, start, end, width, strand, archaic_genotype_summ_Altai_Denisovan, archaic_genotype_summ_Altai_Neanderthal, archaic_genotype_summ_Vindija_Neanderthal, archaic_genotype_summ_Chagyrskaya_Neanderthal, `1KGP_phase3_summ_AF`, `1KGP_phase3_summ_AC`, `1KGP_phase3_summ_AN`, `gnomAD_v3.1.2_summ_AF`, `gnomAD_v3.1.2_summ_AC`, `gnomAD_v3.1.2_summ_AN`, MPRA_VariantID, MPRA_K562_summ_Summary, MPRA_Jurkat_summ_Summary, nearest_gene_summ_Symbol, ABC_summ_Symbol, Roadmap_Epigenomics_summ_Symbol, EpiMap_summ_Symbol, BBJ_finemap_VariantID, BBJ_finemap_orig_method, BBJ_finemap_orig_trait, BBJ_finemap_orig_pip)

# save tables
write_tsv(mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1.tsv.gz"))
write_tsv(mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1.tsv.gz"))

write_tsv(mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_abbrev, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_ENCODE_cCREs_UKBB_finemap_pip0.1_abbrev.tsv.gz"))
write_tsv(mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_abbrev, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_ENCODE_cCREs_BBJ_finemap_pip0.1_abbrev.tsv.gz"))

write_tsv(mpra_emVars_variants_UKBB_finemap_pip0.1, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_UKBB_finemap_pip0.1.tsv.gz"))
write_tsv(mpra_emVars_variants_BBJ_finemap_pip0.1, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_BBJ_finemap_pip0.1.tsv.gz"))

write_tsv(mpra_emVars_variants_UKBB_finemap_pip0.1_abbrev, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_UKBB_finemap_pip0.1_abbrev.tsv.gz"))
write_tsv(mpra_emVars_variants_BBJ_finemap_pip0.1_abbrev, gzfile("../../results/3-overlap_emVars_GWAS_cCREs/mpra_emVars_variants_BBJ_finemap_pip0.1_abbrev.tsv.gz"))
