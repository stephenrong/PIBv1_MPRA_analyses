#!/bin/R

library(data.table)
library(tidyverse)

S_Func_UKB_BBJ_Variants <- fread("../../PIBv1_MPRA_final/results/4-table_corehaps_finemapped_UKBB_BBJ_overlap/Supplementary_Table_SX-Corehaps_fine_mapped_in_UKBB_BBJ_in_final_dataset.txt")

names(S_Func_UKB_BBJ_Variants) <- gsub("Alleles ", "", names(S_Func_UKB_BBJ_Variants))
names(S_Func_UKB_BBJ_Variants) <- gsub("FineMapped ", "Fine mapped ", names(S_Func_UKB_BBJ_Variants))
names(S_Func_UKB_BBJ_Variants)[1:5] <- c("Variant ID", "Sprime allele", "Neanderthal class", "Denisovan class", "Ambiguous class")

S_Func_UKB_BBJ_Variants <- S_Func_UKB_BBJ_Variants %>% 
	separate(col=`Variant ID`, into=c("chrom", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
	mutate(pos = as.numeric(pos)) %>% arrange(desc(`Fine mapped source`), chrom, pos, ref, alt) %>% 
	dplyr::select(-c(chrom, pos, ref, alt)) %>% 
	dplyr::select(`Variant ID`, `Sprime allele`, `Neanderthal class`,	`Denisovan class`, `Ambiguous class`, `Fine mapped source`, `Fine mapped trait`, `Fine mapped description`, `Fine mapped method`, `Fine mapped region`, `Fine mapped csid`, `Fine mapped pip`)

fwrite(S_Func_UKB_BBJ_Variants, file="../submission_files/S_Func_UKB_BBJ_Variants.txt", sep="\t")
