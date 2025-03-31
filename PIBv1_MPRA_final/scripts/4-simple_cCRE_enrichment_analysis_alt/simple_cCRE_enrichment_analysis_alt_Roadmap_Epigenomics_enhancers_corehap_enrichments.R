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

# compute Roadmap_Epigenomics enrichments
Roadmap_Epigenomics_list <- c("Thymus", "Sm_Muscle", "Other", "Myosat", "Neurosph", "iPSC", "Mesench", "Muscle", "Heart", "HSC_and_B-cell", "IMR90", "Epithelial", "ES-deriv", "ESC", "Digestive", "Adipose", "Blood_and_T-cell", "Brain")

# control for DHS background
introgressed_variants <- introgressed_variants %>% 
	filter(!is.na(Roadmap_Epigenomics_enhancers_summ_Summary))

# loop over pop
Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb_list <- list()
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
	Roadmap_Epigenomics_enhancers_corehap_enrichments_class_tb_list <- list()
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
		Roadmap_Epigenomics_enhancers_corehap_enrichments_component_tb_list <- list()
		for (i in Roadmap_Epigenomics_list) {
			print(i)

			col <- paste0("Roadmap_Epigenomics_enhancers_summ_Group-", i)
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp <- introgressed_variants_class_temp
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp$Roadmap_Epigenomics_overlap <- !is.na(introgressed_variants_class_temp[[col]])
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp <- Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp %>% dplyr::select(Is_corehap_variant, Roadmap_Epigenomics_overlap)
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_tab <- table(Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp)
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_fisher <- fisher.test(Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_tab)
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_fisher_tb <- as_tibble(t(c(
				"Fold change"=Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_fisher$estimate[["odds ratio"]], 
				"Number of overlap"=sum(Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp$Roadmap_Epigenomics_overlap), 
				"P-value"=Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_fisher$p.value
			)))
			Roadmap_Epigenomics_enhancers_corehap_enrichments_component_tb_list[[i]] <- Roadmap_Epigenomics_enhancers_corehap_enrichments_component_temp_fisher_tb
		}
		Roadmap_Epigenomics_enhancers_corehap_enrichments_component_tb <- Roadmap_Epigenomics_enhancers_corehap_enrichments_component_tb_list %>% 
			bind_rows(.id = "Group")
		Roadmap_Epigenomics_enhancers_corehap_enrichments_class_tb_list[[j]] <- Roadmap_Epigenomics_enhancers_corehap_enrichments_component_tb
	}
	Roadmap_Epigenomics_enhancers_corehap_enrichments_class_tb <- Roadmap_Epigenomics_enhancers_corehap_enrichments_class_tb_list %>% 
		bind_rows(.id = "Class") %>% 
		mutate(`Class` = factor(`Class`, levels=c("All", "Neanderthal", "Denisovan", "Ambiguous")))
	Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb_list[[k]] <- Roadmap_Epigenomics_enhancers_corehap_enrichments_class_tb
}

Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb <- Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb_list %>% 
	bind_rows(.id = "Pop") %>% 
	mutate(`Pop`=map_pop_to_col_name[`Pop`]) %>% 
	mutate(`Group`=gsub("_", " ", `Group`))

Roadmap_Epigenomics_enhancers_corehap_enrichments_allpop_tb <- Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb %>% 
	filter(Pop == "All") %>% 
	mutate(`FDR-value` = p.adjust(`P-value`, method="BH")) %>%  
	mutate(`Signif` = (`FDR-value` < 0.1))

Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb <- Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb %>% 
	filter(Pop != "All") %>% 
	filter(`Class` == "All") %>% 
	mutate(`FDR-value` = p.adjust(`P-value`, method="BH")) %>%  
	mutate(`Signif` = (`FDR-value` < 0.1))

write_tsv(Roadmap_Epigenomics_enhancers_corehap_enrichments_allpop_tb, gzfile("../../results/4-simple_cCRE_enrichment_analysis_alt/Roadmap_Epigenomics_enhancers_corehap_enrichments_allpop_class_component_tb.txt.gz"))
write_tsv(Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb, gzfile("../../results/4-simple_cCRE_enrichment_analysis_alt/Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_class_component_tb.txt.gz"))

# visualize Roadmap_Epigenomics enrichments
ggplot(Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_tb, 
	aes(x=`Pop`, 
		y=reorder(Group, `Fold change`), 
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
	ggtitle("Roadmap Epigenomics enhancers") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=1,
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/4-simple_cCRE_enrichment_analysis_alt/Roadmap_Epigenomics_enhancers_corehap_enrichments_pop_class_component_tb.pdf", scale=1)

# visualize Roadmap_Epigenomics enrichments
ggplot(Roadmap_Epigenomics_enhancers_corehap_enrichments_allpop_tb, 
	aes(x=`Class`, 
		y=reorder(Group, `Fold change`), 
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
	ggtitle("Roadmap Epigenomics enhancers") + 
	theme(
		axis.text.x=element_text(angle=45, hjust=1), 
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		aspect.ratio=16/4, 
		plot.title = element_text(hjust = 0.5)
	)
ggsave("../../results/4-simple_cCRE_enrichment_analysis_alt/Roadmap_Epigenomics_enhancers_corehap_enrichments_allpop_class_component_tb.pdf", scale=0.8)
