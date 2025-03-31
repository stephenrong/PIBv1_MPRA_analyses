#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(scales)

# for aesthetics
library(RColorBrewer)
library(hrbrthemes)
library(viridis)
library(ggrepel)

# load variant annotations
corehaps_variants_all_annotations_summ <- as_tibble(fread("../../results/3-merge_all_variant_annotations/corehaps_variants_all_annotations_summ.txt.gz"))

corehaps_variants_all_annotations_summ_Ensembl_VEP_full <- 
	corehaps_variants_all_annotations_summ %>% 
	dplyr::select(Functional_Variant_Overall) %>% table() %>% enframe() %>% 
	mutate(percent = 100*value/sum(value)) %>% 
	mutate(name = factor(name, levels=(c(
		"exonic/splicing only", 
		"exonic/splicing+cCRE", 
		"cCRE-only", 
		"other"
	))))

palette <- c("#31a354", "#f2a93b", "#3182bd", "#e7e7e7")
ggplot(corehaps_variants_all_annotations_summ_Ensembl_VEP_full) + 
	aes(0.4, value, fill = name, label = name) + 
	geom_bar(
		position = "stack", 
		stat = "identity", 
		color = "white", 
		linewidth = 0.2
	) + 
	geom_text(
		aes(0.4, value, color = name, 
			label = name
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 1,
		vjust = -0.1,
		size = 2.5,
		angle = 0
	) + 
	geom_text(
		aes(0.4, value, color = name, 
			label = format(value, nsmall=1, big.mark=",")
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.5, 
		vjust = 1.1,
		size = 2,
		angle = 0
	) + 
	scale_fill_manual(name = "", values = palette) + 
	scale_color_manual(name = "", values = palette) + 
	theme_ipsum(base_family = "Arial") + 
	scale_x_discrete(labels = NULL, breaks = NULL) + 
	scale_y_continuous(breaks = NULL) + 
	labs(x = "", y = "") + 
	theme(
		aspect.ratio = 2, 
		legend.position = "none", 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", color = NA),  
		plot.background = element_rect(fill = "transparent", color = NA)
	)
ggsave("../../results/3-visualize_summary_of_functional_consequences/corehaps_variants_all_annotations_summ_Ensembl_VEP_full.pdf", device = cairo_pdf, scale=0.4)


# plot Ensembl VEP other
corehaps_variants_all_annotations_summ_Ensembl_VEP_other <- 
	corehaps_variants_all_annotations_summ %>% 
	filter(!is.na(Functional_Variant_Ensembl_VEP)) %>% 
	dplyr::select(Functional_Variant_Ensembl_VEP) %>% table() %>% enframe() %>% 
	mutate(percent = 100*value/sum(value)) %>% 
	mutate(name = factor(name, levels=(c(
		"3' UTR", 
		"5' UTR", 
		"synonymous", 
		"missense", 
		"splice other", 
		"splice site", 
		"start/stop"
	))))

palette <- (c("#edf8e9", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32"))
ggplot(corehaps_variants_all_annotations_summ_Ensembl_VEP_other) + 
	aes(0.4, value, fill = name, label = name) + 
	geom_bar(
		position = "stack", 
		stat = "identity", 
		color = "white", 
		linewidth = 0.2
	) + 
	coord_flip(clip = "off") + 
	geom_text(
		aes(1, value, color = name, 
			label = name
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.1, 
		vjust = -0.1,
		size = 3,
		angle = 45
	) + 
	geom_text(
		aes(0.4, value, color = name, 
			label = format(value, nsmall=1, big.mark=",")
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.5, 
		vjust = 0.5,
		size = 2,
		angle = 90
	) + 
	scale_fill_manual(name = "", values = palette) + 
	scale_color_manual(name = "", values = palette) + 
	theme_ipsum(base_family = "Arial") + 
	scale_x_discrete(labels = NULL, breaks = NULL) + 
	scale_y_continuous(breaks = NULL) + 
	labs(x = "", y = "") + 
	theme(
		aspect.ratio = 0.33, 
		legend.position = "none", 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", color = NA),  
		plot.background = element_rect(fill = "transparent", color = NA)
	)
ggsave("../../results/3-visualize_summary_of_functional_consequences/corehaps_variants_all_annotations_summ_Ensembl_VEP_other.pdf", device = cairo_pdf, scale = 0.6)


# plot Ensembl VEP other
corehaps_variants_all_annotations_summ_ENCODE_cCREs_other <- 
	corehaps_variants_all_annotations_summ %>% 
	filter(!is.na(Functional_Variant_ENCODE_cCREs)) %>% 
	dplyr::select(Functional_Variant_ENCODE_cCREs) %>% table() %>% enframe() %>% 
	mutate(percent = 100*value/sum(value)) %>% 
	mutate(name = factor(name, levels=(c(
		"TF",
		"CA", 
		"dELS", 
		"pELS",
		"PLS"
	))))

palette <- (c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"))
ggplot(corehaps_variants_all_annotations_summ_ENCODE_cCREs_other) + 
	aes(0.4, value, fill = name, label = name) + 
	geom_bar(
		position = "stack", 
		stat = "identity", 
		color = "white", 
		linewidth = 0.2
	) + 
	coord_flip(clip = "off") + 
	# coord_polar("y", start=pi/2, clip="off") + 
	geom_text(
		aes(1, value, color = name, 
			label = name
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.1, 
		vjust = -0.1,
		size = 3,
		angle = 45
	) + 
	geom_text(
		aes(0.4, value, color = name, 
			label = format(value, nsmall=1, big.mark=",")
		), 
		color = "black", 
		position = position_stack(vjust = 0.5),
		hjust = 0.5, 
		vjust = 0.5,
		size = 2,
		angle = 90
	) + 
	scale_fill_manual(name = "", values = palette) + 
	scale_color_manual(name = "", values = palette) + 
	theme_ipsum(base_family = "Arial") + 
	scale_x_discrete(labels = NULL, breaks = NULL) + 
	scale_y_continuous(breaks = NULL) + 
	labs(x = "", y = "") + 
	theme(
		aspect.ratio = 0.33, 
		legend.position = "none", 
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", color = NA),  
		plot.background = element_rect(fill = "transparent", color = NA)
	)
ggsave("../../results/3-visualize_summary_of_functional_consequences/corehaps_variants_all_annotations_summ_ENCODE_cCREs_other.pdf", device = cairo_pdf, scale = 0.6)
