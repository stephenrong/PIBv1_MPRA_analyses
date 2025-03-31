#!/bin/R

# load packages
library(tidyverse)
library(data.table)
library(eulerr)
library(ggpubr)
library(scales)

# load bio packages
library(plyranges)
library(rtracklayer)

# load files
corehaps_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/corehaps_variants_all_annotations_summ.txt.gz"))

corehaps_variants_finemapped_UKBB <- corehaps_variants %>% 
	filter(!is.na(UKBB_finemap_summ_Summary)) %>% 
	dplyr::select(VariantID, Alleles_SprimeAllele, Alleles_NeandertalClass, Alleles_DenisovanClass, Alleles_AmbiguousClass, UKBB_finemap_summ_Description) %>% 
	separate_longer_delim(cols=c(UKBB_finemap_summ_Description), delim=";") %>% 
	mutate(FineMapped_source="UKBB") %>% 
	separate(col=UKBB_finemap_summ_Description, into=c("FineMapped_trait", "FineMapped_method", "FineMapped_region", "FineMapped_csid", "FineMapped_pip"), sep=",") %>% 
	mutate(FineMapped_pip=as.numeric(gsub("PIP=", "", FineMapped_pip)))

corehaps_variants_finemapped_BBJ <- corehaps_variants %>% 
	filter(!is.na(BBJ_finemap_summ_Summary)) %>% 
	dplyr::select(VariantID, Alleles_SprimeAllele, Alleles_NeandertalClass, Alleles_DenisovanClass, Alleles_AmbiguousClass, BBJ_finemap_summ_Description) %>% 
	separate_longer_delim(cols=c(BBJ_finemap_summ_Description), delim=";") %>% 
	mutate(FineMapped_source="BBJ") %>% 
	separate(col=BBJ_finemap_summ_Description, into=c("FineMapped_trait", "FineMapped_method", "FineMapped_region", "FineMapped_csid", "FineMapped_pip"), sep=",") %>% 
	mutate(FineMapped_pip=as.numeric(gsub("PIP=", "", FineMapped_pip)))

corehaps_variants_finemapped <- bind_rows(corehaps_variants_finemapped_UKBB, corehaps_variants_finemapped_BBJ) %>% 
	arrange(-FineMapped_pip) %>% 
	mutate(FineMapped_duplicated = duplicated(VariantID))

trait_map <- read_tsv("../../../../Datasets/fine_mapped_regions_gwas/Kanai_et_2021_fine_mapped_ukbb/data_download/release1.1/UKBB_94traits_release1.traits")
trait_map <- trait_map %>% dplyr::select(trait, description)
names(trait_map) <- c("FineMapped_trait", "FineMapped_description")

corehaps_variants_finemapped <- left_join(corehaps_variants_finemapped, trait_map)

corehaps_variants_finemapped <- left_join(corehaps_variants_finemapped, trait_map) %>% 
	mutate(FineMapped_description = ifelse(is.na(FineMapped_description), FineMapped_trait, FineMapped_description))

corehaps_variants_finemapped_print <- corehaps_variants_finemapped
names(corehaps_variants_finemapped_print) <- gsub("_", " ", names(corehaps_variants_finemapped_print))

write_tsv(corehaps_variants_finemapped, gzfile("../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/corehaps_variants_finemapped.txt.gz"))
write_tsv(corehaps_variants_finemapped_print, "../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/Supplementary_Table_SX-Corehaps_fine_mapped_in_UKBB_BBJ_in_final_dataset.txt")

corehaps_variants_finemapped_summ <- corehaps_variants_finemapped %>% 
	filter(!FineMapped_duplicated) %>% 
	count(FineMapped_description)
write_tsv(corehaps_variants_finemapped_summ, gzfile("../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/corehaps_variants_finemapped_summ.txt.gz"))

ggplot() + 
  geom_bar(aes(reorder(FineMapped_description, n), n), data=corehaps_variants_finemapped_summ, stat="identity", fill="#4a7cb3") + 
  geom_text(aes(reorder(FineMapped_description, n), n+0.2, label=n), data=corehaps_variants_finemapped_summ, size=2, color="black") + 
  theme_classic() + theme(aspect.ratio = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("#d3d3d3", "#da827f"), name="FDR<0.05") + 
  ylab("High frequency core haplotype SNVs + PIP>0.1") + 
  xlab("UKBB/BBJ trait")
ggsave("../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/corehaps_variants_finemapped_summ.pdf")


# load files
introgressed_variants <- as_tibble(fread("../../results/3-merge_all_variant_annotations/introgressed_variants_all_annotations_summ.txt.gz"))

introgressed_variants_finemapped_UKBB <- introgressed_variants %>% 
	filter(!is.na(UKBB_finemap_summ_Summary)) %>% 
	dplyr::select(VariantID, Alleles_SprimeAllele, Alleles_NeandertalClass, Alleles_DenisovanClass, Alleles_AmbiguousClass, UKBB_finemap_summ_Description) %>% 
	separate_longer_delim(cols=c(UKBB_finemap_summ_Description), delim=";") %>% 
	mutate(FineMapped_source="UKBB") %>% 
	separate(col=UKBB_finemap_summ_Description, into=c("FineMapped_trait", "FineMapped_method", "FineMapped_region", "FineMapped_csid", "FineMapped_pip"), sep=",") %>% 
	mutate(FineMapped_pip=as.numeric(gsub("PIP=", "", FineMapped_pip)))

introgressed_variants_finemapped_BBJ <- introgressed_variants %>% 
	filter(!is.na(BBJ_finemap_summ_Summary)) %>% 
	dplyr::select(VariantID, Alleles_SprimeAllele, Alleles_NeandertalClass, Alleles_DenisovanClass, Alleles_AmbiguousClass, BBJ_finemap_summ_Description) %>% 
	separate_longer_delim(cols=c(BBJ_finemap_summ_Description), delim=";") %>% 
	mutate(FineMapped_source="BBJ") %>% 
	separate(col=BBJ_finemap_summ_Description, into=c("FineMapped_trait", "FineMapped_method", "FineMapped_region", "FineMapped_csid", "FineMapped_pip"), sep=",") %>% 
	mutate(FineMapped_pip=as.numeric(gsub("PIP=", "", FineMapped_pip)))

introgressed_variants_finemapped <- bind_rows(introgressed_variants_finemapped_UKBB, introgressed_variants_finemapped_BBJ) %>% 
	arrange(-FineMapped_pip) %>% 
	mutate(FineMapped_duplicated = duplicated(VariantID))

trait_map <- read_tsv("../../../../Datasets/fine_mapped_regions_gwas/Kanai_et_2021_fine_mapped_ukbb/data_download/release1.1/UKBB_94traits_release1.traits")
trait_map <- trait_map %>% dplyr::select(trait, description)
names(trait_map) <- c("FineMapped_trait", "FineMapped_description")

introgressed_variants_finemapped <- left_join(introgressed_variants_finemapped, trait_map) %>% 
	mutate(FineMapped_description = ifelse(is.na(FineMapped_description), FineMapped_trait, FineMapped_description))

introgressed_variants_finemapped_print <- introgressed_variants_finemapped
names(introgressed_variants_finemapped_print) <- gsub("_", " ", names(introgressed_variants_finemapped_print))

write_tsv(introgressed_variants_finemapped, gzfile("../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/introgressed_variants_finemapped.txt.gz"))
write_tsv(introgressed_variants_finemapped_print, "../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/Supplementary_Table_SX-Introgressed_fine_mapped_in_UKBB_BBJ_in_final_dataset.txt")

introgressed_variants_finemapped_summ <- introgressed_variants_finemapped %>% 
	filter(!FineMapped_duplicated) %>% 
	count(FineMapped_description)
write_tsv(introgressed_variants_finemapped_summ, gzfile("../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/introgressed_variants_finemapped_summ.txt.gz"))

ggplot() + 
  geom_bar(aes(reorder(FineMapped_description, n), n), data=introgressed_variants_finemapped_summ, stat="identity", fill="#4a7cb3") + 
  geom_text(aes(reorder(FineMapped_description, n), n+0.2, label=n), data=introgressed_variants_finemapped_summ, size=2, color="black") + 
  theme_classic() + theme(aspect.ratio = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(values = c("#d3d3d3", "#da827f"), name="FDR<0.05") + 
  ylab("Archaic introgressed SNVs + PIP>0.1") + 
  xlab("UKBB/BBJ trait")
ggsave("../../results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/introgressed_variants_finemapped_summ.pdf", scale=1.5)
