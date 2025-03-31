#!/bin/R

library(data.table)
library(tidyverse)

S_Func_Functional_Annotations_Variants <- fread("../../PIBv1_MPRA_final/results/3-merge_all_variant_annotations/introgressed_variants_all_annotations_summ.txt.gz")
S_Func_Functional_Annotations_Variants <- S_Func_Functional_Annotations_Variants %>% 
	dplyr::select(`seqnames`, `start`, `end`, `width`, `strand`, `VariantID`, `VariantCHROM`, `VariantPOS`, `VariantREF`, `VariantALT`, 
		starts_with("Alleles_"), starts_with("Introgr_TractID_"), starts_with("Introgr_CorehapID_"), starts_with("CorehapID_"), 
		`Flag_variant_failed_liftOver_hg38`, 
		`Ensembl_VEP_summ_Gene`, `Ensembl_VEP_summ_Category`, `Ensembl_VEP_summ_Summary`, `Ensembl_VEP_summ_Symbol`, `SpliceAI_raw_summ_Gene`, `SpliceAI_raw_summ_Description`, `SpliceAI_raw_summ_Score`, `SpliceAI_raw_summ_Summary`, `SpliceAI_raw_summ_Symbol`, `ENCODE_cCREs_summ_Description`, `ENCODE_cCREs_summ_Summary`, 
		`DHS_index_vocabulary_summ_Summary`, `DHS_index_vocabulary_summ_Component`, `Roadmap_Epigenomics_merged_summ_Summary`, `Roadmap_Epigenomics_merged_summ_Group`, `Roadmap_Epigenomics_enhancers_summ_Summary`, `Roadmap_Epigenomics_enhancers_summ_Group`, `Roadmap_Epigenomics_promoters_summ_Summary`, `Roadmap_Epigenomics_promoters_summ_Group`, `Roadmap_Epigenomics_dyadic_summ_Summary`, `Roadmap_Epigenomics_dyadic_summ_Group`, 
	)

names(S_Func_Functional_Annotations_Variants) <- gsub("^Alleles_", "Archaic_", names(S_Func_Functional_Annotations_Variants))
names(S_Func_Functional_Annotations_Variants) <- gsub("^CorehapID_", "HighFreq_CorehapID_", names(S_Func_Functional_Annotations_Variants))

names(S_Func_Functional_Annotations_Variants) <- gsub("^Variant", "Variant_", names(S_Func_Functional_Annotations_Variants))
names(S_Func_Functional_Annotations_Variants) <- gsub("_summ_", "", names(S_Func_Functional_Annotations_Variants))
names(S_Func_Functional_Annotations_Variants) <- gsub("_orig_", "", names(S_Func_Functional_Annotations_Variants))

fwrite(S_Func_Functional_Annotations_Variants, file="../submission_files/S_Func_Functional_Annotations_Variants.txt", sep="\t")
