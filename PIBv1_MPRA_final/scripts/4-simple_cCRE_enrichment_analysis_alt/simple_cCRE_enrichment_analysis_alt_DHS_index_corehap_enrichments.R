#!/bin/R

# load packages
library(tidyverse)
library(data.table)
library(wrapr)

# load bio packages
library(plyranges)

# load files
introgressed_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/introgressed_variants_all_annotations_summ.txt.gz"))

map_pop_to_col_name <- c(
	"Ata"="Ata", 
	"Baining-Kagat"="BainingKagat", 
	"Baining-Mali"="BainingMali", 
	"Bellona"=NA, 
	"Bellona,Rennell"="BellonaRennell", 
	"Goroka"=NA, 
	"Goroka,Sepik"="GorokaSepik", 
	"Kove"="Kove", 
	"Lavongai-Mussau"="LavongaiMussau", 
	"Malaita"="Malaita", 
	"Mamusi"="Mamusi", 
	"Melamela"="Melamela", 
	"Nailik-Notsi-Tigak"="NailikNotsiTigak", 
	"Nakanai"=NA, 
	"Nakanai,Mangseng"="NakanaiMangseng", 
	"Nasioi"="Nasioi", 
	"Santa-Cruz"="SantaCruz", 
	"Saposa"="Saposa", 
	"Sepik"=NA, 
	"Tikopia"="Tikopia", 
	"Vella-Lavella"="VellaLavella"
)

map_pop_to_col_name_short <- map_pop_to_col_name[which(!is.na(map_pop_to_col_name))]

map_pop_to_col_name <- c("All"="All", map_pop_to_col_name)
map_pop_to_col_name <- map_pop_to_col_name[which(!is.na(map_pop_to_col_name))]
map_pop_to_col_name <- setNames(names(map_pop_to_col_name), map_pop_to_col_name)

# filter to those overlapping corehaps
TractID_pop_list <- paste0("Introgr_CorehapID_", map_pop_to_col_name_short)

introgressed_variants <- introgressed_variants %>% 
	mutate(Introgr_Corehaps_count = rowSums(!is.na(across(all_of(TractID_pop_list))))) %>% 
	filter(Introgr_Corehaps_count > 0) %>% 
	dplyr::select(-Introgr_Corehaps_count)

# compute DHS enrichments
DHS_index_list <- c("Lymphoid", "Organ_development_and_renal", "Stromal_B", "Placental_and_trophoblast", "Pulmonary_development", "Cancer_and_epithelial", "Primitive_and_embryonic", "Tissue_invariant", "Cardiac", "Musculoskeletal", "Neural", "Myeloid_and_erythroid", "Renal_and_cancer", "Vascular_and_endothelial", "Digestive", "Stromal_A")

# control for DHS background
introgressed_variants <- introgressed_variants %>% 
	filter(!is.na(DHS_index_vocabulary_summ_Summary))

# loop over pop
DHS_index_corehap_enrichment_pop_tb_list <- list()
for (k in c(names(map_pop_to_col_name))) {
	print(k)

	# get filtered pop tables
	if (k == "All") {
		introgressed_variants_pop_temp <- introgressed_variants %>% 
			# mutate(Is_corehap_variant = (Tracts_count>0))
			mutate(Is_corehap_variant = (Corehaps_count>0))
	} else {
		introgressed_variants_pop_temp <- introgressed_variants %>% 
			filter(!is.na(!!as.name(paste0("Introgr_TractID_", k)))) %>% 
			# mutate(Is_corehap_variant = !is.na(!!as.name(paste0("TractID_", k))))
			mutate(Is_corehap_variant = !is.na(!!as.name(paste0("CorehapID_", k))))
	}

	# loop over class
	DHS_index_corehap_enrichment_class_tb_list <- list()
	for (j in c("All", "Neanderthal", "Denisovan", "Ambiguous")) {
		print(j)

		# get filtered class tables
		if (j == "All") {
			introgressed_variants_class_temp <- introgressed_variants_pop_temp
		} else if (j == "Neanderthal") {
			introgressed_variants_class_temp <- introgressed_variants_pop_temp %>% 
				filter(Alleles_NeandertalClass)
		} else if (j == "Denisovan") {
			introgressed_variants_class_temp <- introgressed_variants_pop_temp %>% 
				filter(Alleles_DenisovanClass)
		} else if (j == "Ambiguous") {
			introgressed_variants_class_temp <- introgressed_variants_pop_temp %>% 
				filter(Alleles_AmbiguousClass)
		}

		# loop over components
		DHS_index_corehap_enrichment_component_tb_list <- list()
		for (i in DHS_index_list) {
			print(i)

			col <- paste0("DHS_index_vocabulary_summ_Component-", i)
			DHS_index_corehap_enrichment_component_temp <- introgressed_variants_class_temp
			DHS_index_corehap_enrichment_component_temp$DHS_index_overlap <- !is.na(introgressed_variants_class_temp[[col]])
			DHS_index_corehap_enrichment_component_temp <- DHS_index_corehap_enrichment_component_temp[c("Is_corehap_variant", "DHS_index_overlap")]
			DHS_index_corehap_enrichment_component_temp_tab <- table(DHS_index_corehap_enrichment_component_temp)
			DHS_index_corehap_enrichment_component_temp_fisher <- fisher.test(DHS_index_corehap_enrichment_component_temp_tab)
			DHS_index_corehap_enrichment_component_temp_fisher_tb <- as_tibble(t(c(
				"Fold change"=DHS_index_corehap_enrichment_component_temp_fisher$estimate[["odds ratio"]], 
				"Number of overlap"=sum(DHS_index_corehap_enrichment_component_temp$DHS_index_overlap), 
				"P-value"=DHS_index_corehap_enrichment_component_temp_fisher$p.value
			)))
			DHS_index_corehap_enrichment_component_tb_list[[i]] <- DHS_index_corehap_enrichment_component_temp_fisher_tb
		}
		DHS_index_corehap_enrichment_component_tb <- DHS_index_corehap_enrichment_component_tb_list %>% 
			bind_rows(.id = "Component")
		DHS_index_corehap_enrichment_class_tb_list[[j]] <- DHS_index_corehap_enrichment_component_tb
	}
	DHS_index_corehap_enrichment_class_tb <- DHS_index_corehap_enrichment_class_tb_list %>% 
		bind_rows(.id = "Class") %>% 
		mutate(`Class` = factor(`Class`, levels=c("All", "Neanderthal", "Denisovan", "Ambiguous")))
	DHS_index_corehap_enrichment_pop_tb_list[[k]] <- DHS_index_corehap_enrichment_class_tb
}

DHS_index_corehap_enrichment_pop_tb <- DHS_index_corehap_enrichment_pop_tb_list %>% 
	bind_rows(.id = "Pop") %>% 
	mutate(`Pop`=map_pop_to_col_name[`Pop`]) %>% 
	mutate(`Component`=gsub("_", " ", `Component`))

DHS_index_corehap_enrichment_allpop_tb <- DHS_index_corehap_enrichment_pop_tb %>% 
	filter(Pop == "All") %>% 
	mutate(`FDR-value` = p.adjust(`P-value`, method="BH")) %>%  
	mutate(`Signif` = (`FDR-value` < 0.1))

DHS_index_corehap_enrichment_pop_tb <- DHS_index_corehap_enrichment_pop_tb %>% 
	filter(Pop != "All") %>% 
	filter(`Class` == "All") %>% 
	mutate(`FDR-value` = p.adjust(`P-value`, method="BH")) %>%  
	mutate(`Signif` = (`FDR-value` < 0.1))

write_tsv(DHS_index_corehap_enrichment_allpop_tb, gzfile("../../results/4-simple_cCRE_enrichment_analysis_alt/DHS_index_corehap_enrichments_allpop_class_component_tb.txt.gz"))
write_tsv(DHS_index_corehap_enrichment_pop_tb, gzfile("../../results/4-simple_cCRE_enrichment_analysis_alt/DHS_index_corehap_enrichments_pop_class_component_tb.txt.gz"))

# visualize DHS enrichments
ggplot(DHS_index_corehap_enrichment_pop_tb, 
	aes(x=`Pop`, 
		y=reorder(Component, `Fold change`), 
		fill=log2(`Fold change`), 
		color=`Signif`,
		size=`Number of overlap`)) + 
	geom_point(
		shape=21,
		stroke=1) + 
	scale_fill_gradient2(low="#D7191D", mid="white", high="#2C7BB6") + 
	scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white")) + 
	labs(y=NULL, x=NULL, color="FDR < 0.1", fill=bquote('Log'[2]~'fold-change'), size="Overlapping variants") + 
	theme_bw() + 
	ggtitle("DHS index components") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=1,  # 16/4, 
		plot.title = element_text(hjust = 0.5)  # , 
	)
ggsave("../../results/4-simple_cCRE_enrichment_analysis_alt/DHS_index_corehap_enrichments_pop_class_component_tb.pdf", scale=1)

# visualize DHS enrichments
ggplot(DHS_index_corehap_enrichment_allpop_tb, 
	aes(x=`Class`, 
		y=reorder(Component, `Fold change`), 
		fill=log2(`Fold change`), 
		color=`Signif`,
		size=`Number of overlap`)) + 
	geom_point(
		shape=21,
		stroke=1) + 
	scale_fill_gradient2(low="#D7191D", mid="white", high="#2C7BB6") + 
	scale_color_manual(breaks=c(TRUE, FALSE), values=c("black", "white")) + 
	labs(y=NULL, x=NULL, color="FDR < 0.1", fill=bquote('Log'[2]~'fold-change'), size="Overlapping variants") + 
	theme_bw() + 
	ggtitle("DHS index components") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=16/4, 
		plot.title = element_text(hjust = 0.5)  # , 
	)
ggsave("../../results/4-simple_cCRE_enrichment_analysis_alt/DHS_index_corehap_enrichments_allpop_class_component_tb.pdf", scale=0.8)
