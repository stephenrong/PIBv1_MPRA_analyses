#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# load summ files
# 	all variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz"))

# 	population allele frequencies
adaptive_variants_1KGP_phase3_AFs_summ <- as_tibble(fread("../../results/2c-annotate_variants_1KGP_AFs/adaptive_variants_1KGP_phase3_AFs_summ.txt.gz"))
adaptive_variants_gnomAD_v3.1.2_AFs_summ <- as_tibble(fread("../../results/2d-annotate_variants_gnomAD_AFs/adaptive_variants_gnomAD_v3.1.2_AFs_summ.txt.gz"))
adaptive_variants_archaic_genotypes_summ <- as_tibble(fread("../../results/2e-annotate_variants_archaic_genotypes/adaptive_variants_archaic_genotypes_summ.txt.gz"))

# 	coding and splicing
adaptive_variants_Ensembl_VEP_summ <- as_tibble(fread("../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants_Ensembl_VEP_summ.txt.gz"))
adaptive_variants_SpliceAI_summ <- as_tibble(fread("../../results/2b-annotate_variants_SpliceAI/adaptive_variants_SpliceAI_summ.txt.gz"))

# 	evolutionary conservation
adaptive_variants_CADD_summ <- as_tibble(fread("../../results/2f-annotate_variants_CADD_scores/adaptive_variants_CADD_summ.txt.gz"))
adaptive_variants_phyloP_summ <- as_tibble(fread("../../results/2g-annotate_variants_phyloP_scores/adaptive_variants_phyloP_summ.txt.gz"))

# 	miscellaneous annotations
adaptive_variants_ClinVar_summ <- as_tibble(fread("../../results/2h-annotate_variants_ClinVar_categories/adaptive_variants_ClinVar_summ.txt.gz"))
adaptive_variants_BBJ_finemap_summ <- as_tibble(fread("../../results/2m-annotate_variants_UKBB_BBJ_finemap/adaptive_variants_BBJ_finemap_summ.txt.gz"))
adaptive_variants_UKBB_finemap_summ <- as_tibble(fread("../../results/2m-annotate_variants_UKBB_BBJ_finemap/adaptive_variants_UKBB_finemap_summ.txt.gz"))

# 	MPRA emVar annotations
adaptive_variants_MPRA_K562_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562_summ.txt.gz"))
adaptive_variants_MPRA_Jurkat_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat_summ.txt.gz"))
adaptive_variants_MPRA_HepG2_summ <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_HepG2_summ.txt.gz"))

# 	ENCODE cCRE overlap
adaptive_variants_ENCODE_cCREs_summ <- as_tibble(fread("../../results/2j-annotate_variants_ENCODE_cCREs/adaptive_variants_ENCODE_cCREs_summ.txt.gz"))

# 	ENCODE cell-line cCRE overlap
adaptive_variants_ENCODE_cCREs_K562_summ <- as_tibble(fread("../../results/2k-annotate_variants_ENCODE_cell_lines/adaptive_variants_ENCODE_cCREs_K562_summ.txt.gz"))
adaptive_variants_ENCODE_cCREs_Jurkat_summ <- as_tibble(fread("../../results/2k-annotate_variants_ENCODE_cell_lines/adaptive_variants_ENCODE_cCREs_Jurkat_summ.txt.gz"))
adaptive_variants_ENCODE_cCREs_HepG2_summ <- as_tibble(fread("../../results/2k-annotate_variants_ENCODE_cell_lines/adaptive_variants_ENCODE_cCREs_HepG2_summ.txt.gz"))

# 	DHS vocabulary
adaptive_variants_DHS_index_vocabulary_summ <- as_tibble(fread("../../results/2r-annotate_variants_DHS_index_vocabulary/adaptive_variants_DHS_index_vocabulary_summ.txt.gz"))

# 	Roadmap Epigenomics groups
adaptive_variants_Roadmap_Epigenomics_merged_summ <- as_tibble(fread("../../results/2q-annotate_variants_Roadmap_Epigenomics_groups/adaptive_variants_Roadmap_Epigenomics_merged_summ.txt.gz"))
adaptive_variants_Roadmap_Epigenomics_enhancers_summ <- as_tibble(fread("../../results/2q-annotate_variants_Roadmap_Epigenomics_groups/adaptive_variants_Roadmap_Epigenomics_enhancers_summ.txt.gz"))
adaptive_variants_Roadmap_Epigenomics_promoters_summ <- as_tibble(fread("../../results/2q-annotate_variants_Roadmap_Epigenomics_groups/adaptive_variants_Roadmap_Epigenomics_promoters_summ.txt.gz"))
adaptive_variants_Roadmap_Epigenomics_dyadic_summ <- as_tibble(fread("../../results/2q-annotate_variants_Roadmap_Epigenomics_groups/adaptive_variants_Roadmap_Epigenomics_dyadic_summ.txt.gz"))

# 	TF motif annotations
adaptive_variants_TF_motifs_hocomoco_summ <- as_tibble(fread("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_hocomoco_summ.txt.gz"))
adaptive_variants_TF_motifs_jaspar_summ <- as_tibble(fread("../../results/2n-annotate_variants_TF_motifs/adaptive_variants_TF_motifs_jaspar_summ.txt.gz"))

# 	DNase TF footprints
adaptive_variants_DNase_footprints_summ <- as_tibble(fread("../../results/2s-annotate_variants_DNase_footprints/adaptive_variants_DNase_footprints_summ.txt.gz"))

# 	ReMap motif binding sites
adaptive_variants_ReMap_TFBS_K562_summ <- as_tibble(fread("../../results/2l-annotate_variants_ReMap_TFBS/adaptive_variants_ReMap_TFBS_K562_summ.txt.gz"))
adaptive_variants_ReMap_TFBS_Jurkat_summ <- as_tibble(fread("../../results/2l-annotate_variants_ReMap_TFBS/adaptive_variants_ReMap_TFBS_Jurkat_summ.txt.gz"))
adaptive_variants_ReMap_TFBS_HepG2_summ <- as_tibble(fread("../../results/2l-annotate_variants_ReMap_TFBS/adaptive_variants_ReMap_TFBS_HepG2_summ.txt.gz"))

# 	Nearest gene
adaptive_variants_nearest_gene_summ <- as_tibble(fread("../../results/2o-annotate_variants_nearest_gene/adaptive_variants_nearest_gene_summ.txt.gz"))

# 	Gene enhancer links
adaptive_variants_ABC_summ <- as_tibble(fread("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_ABC_summ.txt.gz"))
adaptive_variants_Roadmap_Epigenomics_summ <- as_tibble(fread("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_Roadmap_Epigenomics_summ.txt.gz"))
adaptive_variants_EpiMap_summ <- as_tibble(fread("../../results/2p-annotate_variants_gene_enhancer_links/adaptive_variants_EpiMap_summ.txt.gz"))

# merge summ files
adaptive_variants_all_annotations_summ <- adaptive_variants_tb %>% 
	left_join(adaptive_variants_archaic_genotypes_summ) %>% 
	left_join(adaptive_variants_1KGP_phase3_AFs_summ) %>% 
	left_join(adaptive_variants_gnomAD_v3.1.2_AFs_summ) %>% 
	left_join(adaptive_variants_Ensembl_VEP_summ) %>% 
	left_join(adaptive_variants_SpliceAI_summ) %>% 
	left_join(adaptive_variants_CADD_summ) %>% 
	left_join(adaptive_variants_phyloP_summ) %>% 
	left_join(adaptive_variants_ClinVar_summ) %>% 
	left_join(adaptive_variants_BBJ_finemap_summ) %>% 
	left_join(adaptive_variants_UKBB_finemap_summ) %>% 
	left_join(adaptive_variants_MPRA_K562_summ) %>% 
	left_join(adaptive_variants_MPRA_Jurkat_summ) %>% 
	left_join(adaptive_variants_ENCODE_cCREs_summ) %>% 
	left_join(adaptive_variants_ENCODE_cCREs_K562_summ) %>% 
	left_join(adaptive_variants_ENCODE_cCREs_Jurkat_summ) %>% 
	left_join(adaptive_variants_TF_motifs_hocomoco_summ) %>% 
	left_join(adaptive_variants_TF_motifs_jaspar_summ) %>% 
	left_join(adaptive_variants_DNase_footprints_summ) %>% 
	left_join(adaptive_variants_ReMap_TFBS_K562_summ) %>% 
	left_join(adaptive_variants_ReMap_TFBS_Jurkat_summ) %>% 
	left_join(adaptive_variants_nearest_gene_summ) %>% 
	left_join(adaptive_variants_DHS_index_vocabulary_summ) %>% 
	left_join(adaptive_variants_Roadmap_Epigenomics_merged_summ) %>% 
	left_join(adaptive_variants_Roadmap_Epigenomics_enhancers_summ) %>% 
	left_join(adaptive_variants_Roadmap_Epigenomics_promoters_summ) %>% 
	left_join(adaptive_variants_Roadmap_Epigenomics_dyadic_summ) %>% 
	left_join(adaptive_variants_ABC_summ) %>% 
	left_join(adaptive_variants_Roadmap_Epigenomics_summ) %>% 
	left_join(adaptive_variants_EpiMap_summ)

# get top level functional annotation categories
Ensembl_VEP_full_map <- c(
	"3_prime_UTR_variant"="exonic/splicing only", 
	"5_prime_UTR_variant"="exonic/splicing only", 
	"coding_sequence_variant"="other",
	"downstream_gene_variant"="other", 
	"intergenic_variant"="other", 
	"intron_variant"="other", 
	"mature_miRNA_variant"="other", 
	"missense_variant"="exonic/splicing only", 
	"non_coding_transcript_exon_variant"="other", 
	"splice_acceptor_variant"="exonic/splicing only", 
	"splice_donor_5th_base_variant"="exonic/splicing only", 
	"splice_donor_region_variant"="exonic/splicing only", 
	"splice_donor_variant"="exonic/splicing only", 
	"splice_polypyrimidine_tract_variant"="exonic/splicing only", 
	"splice_region_variant"="exonic/splicing only", 
	"start_gained"="exonic/splicing only", 
	"start_lost"="exonic/splicing only", 
	"stop_gained"="exonic/splicing only", 
	"stop_lost"="exonic/splicing only", 
	"synonymous_variant"="exonic/splicing only", 
	"stop_retained_variant"="other",
	"upstream_gene_variant"="other"
)

adaptive_variants_all_annotations_summ <- 
	adaptive_variants_all_annotations_summ %>% 
	mutate(Functional_Variant_Overall = Ensembl_VEP_full_map[Ensembl_VEP_summ_Category]) %>% 
	mutate(Functional_Variant_Overall = ifelse((!is.na(ENCODE_cCREs_summ_Summary) & (Functional_Variant_Overall == "other")), "cCRE-only", Functional_Variant_Overall)) %>% 
	mutate(Functional_Variant_Overall = ifelse((!is.na(ENCODE_cCREs_summ_Summary) & (Functional_Variant_Overall == "exonic/splicing only")), "exonic/splicing+cCRE", Functional_Variant_Overall))

Ensembl_VEP_other_map <- c(
	"3_prime_UTR_variant"="3' UTR", 
	"5_prime_UTR_variant"="5' UTR", 
	"coding_sequence_variant"=NA,
	"downstream_gene_variant"=NA, 
	"intergenic_variant"=NA, 
	"intron_variant"=NA, 
	"mature_miRNA_variant"=NA, 
	"missense_variant"="missense", 
	"non_coding_transcript_exon_variant"=NA, 
	"splice_acceptor_variant"="splice site", 
	"splice_donor_5th_base_variant"="splice other", 
	"splice_donor_region_variant"="splice other", 
	"splice_donor_variant"="splice site", 
	"splice_polypyrimidine_tract_variant"="splice other", 
	"splice_region_variant"="splice other", 
	"start_gained"="start/stop", 
	"start_lost"="start/stop",
	"stop_gained"="start/stop",
	"stop_lost"="start/stop",
	"synonymous_variant"="synonymous", 
	"stop_retained_variant"=NA,
	"upstream_gene_variant"=NA
)

adaptive_variants_all_annotations_summ <- 
	adaptive_variants_all_annotations_summ %>% 
	mutate(Functional_Variant_Ensembl_VEP = Ensembl_VEP_other_map[Ensembl_VEP_summ_Category])

ENCODE_cCREs_other_map <- c(
	"CA"="CA",
	"CA-CTCF"="CA",
	"CA-H3K4me3"="CA",
	"CA-TF"="CA",
	"CA-TF,CA"="CA",
	"CA,CA"="CA",
	"dELS"="dELS",
	"dELS,CA"="dELS",
	"dELS,CA-CTCF"="dELS",
	"dELS,dELS"="dELS",
	"pELS"="pELS",
	"PLS"="PLS",
	"TF"="TF"
)

adaptive_variants_all_annotations_summ <- 
	adaptive_variants_all_annotations_summ %>% 
	mutate(Functional_Variant_ENCODE_cCREs = ENCODE_cCREs_other_map[ENCODE_cCREs_summ_Summary])

# save merge summ
write_tsv(adaptive_variants_all_annotations_summ, gzfile("../../results/3-merge_all_variant_annotations/adaptive_variants_all_annotations_summ.txt.gz"))

mpra_ccre_variants_all_annotations_summ <- adaptive_variants_all_annotations_summ %>% filter(PIB_pilot_mpra_ccre_variant)
write_tsv(mpra_ccre_variants_all_annotations_summ, gzfile("../../results/3-merge_all_variant_annotations/mpra_ccre_variants_all_annotations_summ.txt.gz"))
