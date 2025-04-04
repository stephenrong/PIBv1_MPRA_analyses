#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load TSV
adaptive_variants_MPRA_Jurkat <- as_tibble(fread("../../PIBv1_MPRA_mpra/Tucci_MPRA_pilot/MPRAmodel_output_Jurkat/results/MPRAmodel_PIBv1_MPRA_pilot_Jurkat_Jurkat_emVAR_glm_20230106.out"))

# add prefix and save names
names(adaptive_variants_MPRA_Jurkat) <- paste("MPRA_Jurkat_orig", names(adaptive_variants_MPRA_Jurkat), sep="_")
adaptive_variants_MPRA_Jurkat_names <- names(adaptive_variants_MPRA_Jurkat)

# add GRange and Variant cols
adaptive_variants_MPRA_Jurkat <- adaptive_variants_MPRA_Jurkat %>% 
	mutate(VariantID = MPRA_Jurkat_orig_SNP) %>% 
	separate(col=VariantID, into=c("VariantCHROM", "VariantPOS", "VariantREF", "VariantALT"), sep="_", remove=FALSE) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, all_of(adaptive_variants_MPRA_Jurkat_names))

# add summary cols
adaptive_variants_MPRA_Jurkat <- adaptive_variants_MPRA_Jurkat %>% 
	mutate(MPRA_Jurkat_summ_Description_activity_log10FDR = pmax(MPRA_Jurkat_orig_A_logPadj_BH, MPRA_Jurkat_orig_B_logPadj_BH)) %>% 
	mutate(MPRA_Jurkat_summ_Description_activity_log2FC = pmax(MPRA_Jurkat_orig_A_log2FC, MPRA_Jurkat_orig_B_log2FC)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log10FDR = as.numeric(MPRA_Jurkat_orig_Skew_logFDR_act)) %>% 
	mutate(MPRA_Jurkat_summ_Description_skew_log2FC = as.numeric(MPRA_Jurkat_orig_Log2Skew)) %>% 
	mutate(MPRA_Jurkat_summ_Description_active = ifelse(
		MPRA_Jurkat_summ_Description_activity_log10FDR > -log10(0.01), "Jurkat_active", NA)  # FDR < 0.01 for max allele
	) %>% 
	mutate(MPRA_Jurkat_summ_Summary = 
		ifelse((!is.na(MPRA_Jurkat_summ_Description_active)) & (MPRA_Jurkat_summ_Description_skew_log10FDR > -log10(0.05)), # FDR < 0.05
		ifelse(MPRA_Jurkat_summ_Description_skew_log2FC>0, "Jurkat_pos_emVar", "Jurkat_neg_emVar"), NA)  # distinguish positive negative
	)

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(adaptive_variants_tb), names(adaptive_variants_MPRA_Jurkat))
adaptive_variants_MPRA_Jurkat <- adaptive_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(adaptive_variants_MPRA_Jurkat %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(adaptive_variants_MPRA_Jurkat, gzfile("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat.txt.gz"))

# save summ variants
adaptive_variants_MPRA_Jurkat_summ <- adaptive_variants_MPRA_Jurkat %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_MPRA_Jurkat_summ, gzfile("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat_summ.txt.gz"))

# visual check
pdf("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat_summ_logFDR_sort.pdf")
plot(sort(adaptive_variants_MPRA_Jurkat$MPRA_Jurkat_summ_Description_skew_log10FDR), ylab="MPRA Jurkat logFDR")
dev.off()

pdf("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat_summ_skew_log2FC_sort.pdf")
plot(sort(adaptive_variants_MPRA_Jurkat$MPRA_Jurkat_summ_Description_skew_log2FC), ylab="MPRA Jurkat logSkew")
dev.off()
