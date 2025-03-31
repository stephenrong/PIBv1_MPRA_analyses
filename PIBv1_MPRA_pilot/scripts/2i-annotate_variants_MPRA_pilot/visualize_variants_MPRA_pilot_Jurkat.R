#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
source("../shared_functions/seqinfo_fix_change.R")

# load motifbreakR input
adaptive_variants_MPRA_Jurkat <- as_tibble(fread("../../results/2i-annotate_variants_MPRA_pilot/adaptive_variants_MPRA_Jurkat.txt.gz"))

# temporary filters
adaptive_variants_MPRA_Jurkat <- adaptive_variants_MPRA_Jurkat %>% 
	filter(VariantID %in% c(
		"12_113346546_C_G", "12_113347017_T_C", "12_113348921_G_T", "12_113348935_G_C",  # ,  # OAS1
		"1_89587575_T_C", "1_89592457_C_T", "1_89592489_T_C", "1_89594982_T_G", "1_89596265_G_A", "1_89599705_A_G",  # GBP2
		"6_138020511_C_T", "6_138037108_G_A", "6_138039111_A_T", "6_138226361_A_C", "6_138248445_G_A", "6_138786583_G_C", #,  # TNFAIP3
		"8_116327463_A_G", "8_116364887_T_C", "8_116393745_G_A", "8_116394033_G_T",  # TRPS1
		"1_65419060_C_T", "1_65431457_G_A", "1_65475897_G_T", "1_65475965_A_G", "1_65475981_T_A", "1_65468183_T_A", "1_65477791_C_T",  #,  # JAK1
		"9_96576740_C_G",  # BARX1
		"2_161933715_T_C", "2_162016242_G_A", "2_162020918_T_C"  # TANK
		)
	)

# visualize
for (i in adaptive_variants_MPRA_Jurkat$VariantID) {
	adaptive_variants_MPRA_Jurkat_temp <- adaptive_variants_MPRA_Jurkat %>% 
		filter(VariantID == i)
	ggplot(adaptive_variants_MPRA_Jurkat_temp) + 
		geom_bar(aes(1, MPRA_Jurkat_orig_A_log2FC), stat="identity", fill="#BDD7E7") + 
		geom_errorbar(aes(1, ymin=MPRA_Jurkat_orig_A_log2FC-2*MPRA_Jurkat_orig_A_log2FC_SE, ymax=MPRA_Jurkat_orig_A_log2FC+2*MPRA_Jurkat_orig_A_log2FC_SE), stat="identity", width=0.2) + 
		geom_bar(aes(2, MPRA_Jurkat_orig_B_log2FC), stat="identity", fill="#2171B5") + 
		geom_errorbar(aes(2, ymin=MPRA_Jurkat_orig_B_log2FC-2*MPRA_Jurkat_orig_B_log2FC_SE, ymax=MPRA_Jurkat_orig_B_log2FC+2*MPRA_Jurkat_orig_B_log2FC_SE), stat="identity", width=0.2) + 
		scale_x_discrete("Allele", labels = c("REF", "ALT"), limits = c(1, 2)) + 
		ylab("MPRA Jurkat activity (log2 FC)") + 
		theme(
			aspect.ratio = 2, 
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_rect(fill = "transparent", colour = NA),  
			plot.background = element_rect(fill = "transparent", colour = NA), 
			axis.line = element_line(colour = "black")
		)
	ggsave(paste0("../../results/2i-annotate_variants_MPRA_pilot/variant_level_visualizations/variant_MPRA_Jurkat_activity-", i, ".pdf"), width=3, height=3)
}