#!/bin/R

library(data.table)
library(tidyverse)

S_Func_MPRA_Variants <- fread("../../PIBv1_MPRA_pilot/results/3-merge_all_variant_annotations/mpra_ccre_variants_all_annotations_summ.txt.gz")
S_Func_MPRA_Variants_K562 <- fread("../../PIBv1_MPRA_pilot/results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_K562.txt.gz")
S_Func_MPRA_Variants_Jurkat <- fread("../../PIBv1_MPRA_pilot/results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat.txt.gz")

S_Func_MPRA_Variants <- S_Func_MPRA_Variants %>% 
	left_join(S_Func_MPRA_Variants_K562) %>% 
	left_join(S_Func_MPRA_Variants_Jurkat)

S_Func_MPRA_Variants <- S_Func_MPRA_Variants %>% 
	dplyr::select(`seqnames`, `start`, `end`, `width`, `strand`, `VariantID`, `VariantCHROM`, `VariantPOS`, `VariantREF`, `VariantALT`, 
		`Pop`, `TractID`, `VariantARCALLELE`, 
		`Ata_TractID`, `Baining_TractID`, `Lavongai_TractID`, `Mamusi_TractID`, `Melamela_TractID`, 
		`Flag_variant_failed_liftOver_hg38`, 
		`Ensembl_VEP_summ_Gene`, `Ensembl_VEP_summ_Category`, `Ensembl_VEP_summ_Summary`, `Ensembl_VEP_summ_Symbol`, `SpliceAI_raw_summ_Gene`, `SpliceAI_raw_summ_Description`, `SpliceAI_raw_summ_Score`, `SpliceAI_raw_summ_Summary`, `SpliceAI_raw_summ_Symbol`, `ENCODE_cCREs_summ_Description`, `ENCODE_cCREs_summ_Summary`, 
		`DHS_index_vocabulary_summ_Summary`, `DHS_index_vocabulary_summ_Component`, `Roadmap_Epigenomics_merged_summ_Summary`, `Roadmap_Epigenomics_merged_summ_Group`, `Roadmap_Epigenomics_enhancers_summ_Summary`, `Roadmap_Epigenomics_enhancers_summ_Group`, `Roadmap_Epigenomics_promoters_summ_Summary`, `Roadmap_Epigenomics_promoters_summ_Group`, `Roadmap_Epigenomics_dyadic_summ_Summary`, `Roadmap_Epigenomics_dyadic_summ_Group`, 
		`MPRA_K562_orig_ID`, `MPRA_K562_orig_comb`, `MPRA_K562_orig_SNP`, `MPRA_K562_orig_chr`, `MPRA_K562_orig_pos`, `MPRA_K562_orig_ref_allele`, `MPRA_K562_orig_alt_allele`, `MPRA_K562_orig_allele`, `MPRA_K562_orig_window`, `MPRA_K562_orig_strand`, `MPRA_K562_orig_haplotype`, `MPRA_K562_orig_A_Ctrl_Mean`, `MPRA_K562_orig_A_Exp_Mean`, `MPRA_K562_orig_A_log2FC`, `MPRA_K562_orig_A_log2FC_SE`, `MPRA_K562_orig_A_logP`, `MPRA_K562_orig_A_logPadj_BH`, `MPRA_K562_orig_A_logPadj_BF`, `MPRA_K562_orig_B_Ctrl_Mean`, `MPRA_K562_orig_B_Exp_Mean`, `MPRA_K562_orig_B_log2FC`, `MPRA_K562_orig_B_log2FC_SE`, `MPRA_K562_orig_B_logP`, `MPRA_K562_orig_B_logPadj_BH`, `MPRA_K562_orig_B_logPadj_BF`, `MPRA_K562_orig_Log2Skew`, `MPRA_K562_orig_Skew_SE`, `MPRA_K562_orig_skewStat`, `MPRA_K562_orig_Skew_logP`, `MPRA_K562_orig_Skew_logFDR`, `MPRA_K562_orig_Skew_logFDR_act`, `MPRA_K562_summ_Summary`, 
		`MPRA_Jurkat_orig_ID`, `MPRA_Jurkat_orig_comb`, `MPRA_Jurkat_orig_SNP`, `MPRA_Jurkat_orig_chr`, `MPRA_Jurkat_orig_pos`, `MPRA_Jurkat_orig_ref_allele`, `MPRA_Jurkat_orig_alt_allele`, `MPRA_Jurkat_orig_allele`, `MPRA_Jurkat_orig_window`, `MPRA_Jurkat_orig_strand`, `MPRA_Jurkat_orig_haplotype`, `MPRA_Jurkat_orig_A_Ctrl_Mean`, `MPRA_Jurkat_orig_A_Exp_Mean`, `MPRA_Jurkat_orig_A_log2FC`, `MPRA_Jurkat_orig_A_log2FC_SE`, `MPRA_Jurkat_orig_A_logP`, `MPRA_Jurkat_orig_A_logPadj_BH`, `MPRA_Jurkat_orig_A_logPadj_BF`, `MPRA_Jurkat_orig_B_Ctrl_Mean`, `MPRA_Jurkat_orig_B_Exp_Mean`, `MPRA_Jurkat_orig_B_log2FC`, `MPRA_Jurkat_orig_B_log2FC_SE`, `MPRA_Jurkat_orig_B_logP`, `MPRA_Jurkat_orig_B_logPadj_BH`, `MPRA_Jurkat_orig_B_logPadj_BF`, `MPRA_Jurkat_orig_Log2Skew`, `MPRA_Jurkat_orig_Skew_SE`, `MPRA_Jurkat_orig_skewStat`, `MPRA_Jurkat_orig_Skew_logP`, `MPRA_Jurkat_orig_Skew_logFDR`, `MPRA_Jurkat_orig_Skew_logFDR_act`, `MPRA_Jurkat_summ_Summary`, 
		`ENCODE_cCREs_K562_summ_Category`, `ENCODE_cCREs_K562_summ_Description`, `ENCODE_cCREs_K562_summ_Summary`, `ENCODE_cCREs_Jurkat_summ_Category`, `ENCODE_cCREs_Jurkat_summ_Description`, `ENCODE_cCREs_Jurkat_summ_Summary`, `TF_motifs_hocomoco_summ_Count`, `TF_motifs_hocomoco_summ_Description`, `TF_motifs_hocomoco_summ_Summary`, `TF_motifs_jaspar_summ_Count`, `TF_motifs_jaspar_summ_Description`, `TF_motifs_jaspar_summ_Summary`, `DNase_footprints_summ_Count`, `DNase_footprints_summ_Description`, `DNase_footprints_summ_Summary`, `ReMap_TFBS_K562_summ_Count`, `ReMap_TFBS_K562_summ_Description`, `ReMap_TFBS_K562_summ_Summary`, `ReMap_TFBS_Jurkat_summ_Count`, `ReMap_TFBS_Jurkat_summ_Description`, `ReMap_TFBS_Jurkat_summ_Summary`, 
		`ABC_summ_Summary`, `ABC_summ_Count`, `ABC_summ_Score`, `ABC_summ_Gene`, `ABC_summ_Symbol`, `ABC_summ_Biosample`, `ABC_summ_Category`, `Roadmap_Epigenomics_summ_Summary`, `Roadmap_Epigenomics_summ_Count`, `Roadmap_Epigenomics_summ_Score`, `Roadmap_Epigenomics_summ_Ensembl_ID`, `Roadmap_Epigenomics_summ_Symbol`, `Roadmap_Epigenomics_summ_Biosample`, `Roadmap_Epigenomics_summ_Category`, `EpiMap_summ_Summary`, `EpiMap_summ_Count`, `EpiMap_summ_Score`, `EpiMap_summ_Ensembl_ID`, `EpiMap_summ_Symbol`, `EpiMap_summ_Group`
	)

S_Func_MPRA_Variants <- S_Func_MPRA_Variants %>% 
	dplyr::rename(
		`MPRA_Pop` = `Pop`, 
		`MRPA_TractID` = `TractID`,
		`MPRA_Archaic_SprimeAllele` = `VariantARCALLELE`
	)

names(S_Func_MPRA_Variants) <- gsub("^Variant", "Variant_", names(S_Func_MPRA_Variants))
names(S_Func_MPRA_Variants) <- gsub("_summ_", "", names(S_Func_MPRA_Variants))
names(S_Func_MPRA_Variants) <- gsub("_orig_", "", names(S_Func_MPRA_Variants))

fwrite(S_Func_MPRA_Variants, file="../submission_files/S_Func_MPRA_Variants.txt", sep="\t")
