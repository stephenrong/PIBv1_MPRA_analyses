#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)

# load vcfs
introgressed_variants_altai_denisovan <- read.vcfR("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-altai_denisovan_masked_norm.vcf.gz")
introgressed_variants_altai_neanderthal <- read.vcfR("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-altai_neanderthal_masked_norm.vcf.gz")
introgressed_variants_vindija_neanderthal <- read.vcfR("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-vindija_neanderthal_masked_norm.vcf.gz")
introgressed_variants_chagyrskaya_neanderthal <- read.vcfR("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants-chagyrskaya_neanderthal_masked_norm.vcf.gz")

introgressed_variants_altai_denisovan <- bind_cols(as_tibble(introgressed_variants_altai_denisovan@fix), as_tibble(introgressed_variants_altai_denisovan@gt))
introgressed_variants_altai_neanderthal <- bind_cols(as_tibble(introgressed_variants_altai_neanderthal@fix), as_tibble(introgressed_variants_altai_neanderthal@gt))
introgressed_variants_vindija_neanderthal <- bind_cols(as_tibble(introgressed_variants_vindija_neanderthal@fix), as_tibble(introgressed_variants_vindija_neanderthal@gt))
introgressed_variants_chagyrskaya_neanderthal <- bind_cols(as_tibble(introgressed_variants_chagyrskaya_neanderthal@fix), as_tibble(introgressed_variants_chagyrskaya_neanderthal@gt))

# process vcfs
introgressed_variants_altai_denisovan <- introgressed_variants_altai_denisovan %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	mutate(archaic_genotype_summ_Altai_Denisovan = gsub(":.*", "", Denisova)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, archaic_genotype_summ_Altai_Denisovan)

introgressed_variants_altai_neanderthal <- introgressed_variants_altai_neanderthal %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	mutate(archaic_genotype_summ_Altai_Neanderthal = gsub(":.*", "", AltaiNeandertal)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, archaic_genotype_summ_Altai_Neanderthal)

introgressed_variants_vindija_neanderthal <- introgressed_variants_vindija_neanderthal %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	mutate(archaic_genotype_summ_Vindija_Neanderthal = gsub(":.*", "", `Vindija33.19`)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, archaic_genotype_summ_Vindija_Neanderthal)

introgressed_variants_chagyrskaya_neanderthal <- introgressed_variants_chagyrskaya_neanderthal %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	mutate(archaic_genotype_summ_Chagyrskaya_Neanderthal = gsub(":.*", "", `Chagyrskaya-Phalanx`)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, archaic_genotype_summ_Chagyrskaya_Neanderthal)

# split into ref and alt
introgressed_variants_altai_denisovan_ref <- introgressed_variants_altai_denisovan %>% filter(!is.na(VariantALT)) %>% 
	mutate(archaic_genotype_summ_Altai_Denisovan = ifelse(archaic_genotype_summ_Altai_Denisovan == "1/0", "0/1", archaic_genotype_summ_Altai_Denisovan))
introgressed_variants_altai_denisovan_alt <- introgressed_variants_altai_denisovan %>% filter(is.na(VariantALT)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantCHROM, VariantPOS, archaic_genotype_summ_Altai_Denisovan)

introgressed_variants_altai_neanderthal_ref <- introgressed_variants_altai_neanderthal %>% filter(!is.na(VariantALT)) %>% 
	mutate(archaic_genotype_summ_Altai_Neanderthal = ifelse(archaic_genotype_summ_Altai_Neanderthal == "1/0", "0/1", archaic_genotype_summ_Altai_Neanderthal))
introgressed_variants_altai_neanderthal_alt <- introgressed_variants_altai_neanderthal %>% filter(is.na(VariantALT)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantCHROM, VariantPOS, archaic_genotype_summ_Altai_Neanderthal)

introgressed_variants_vindija_neanderthal_ref <- introgressed_variants_vindija_neanderthal %>% filter(!is.na(VariantALT)) %>% 
	mutate(archaic_genotype_summ_Vindija_Neanderthal = ifelse(archaic_genotype_summ_Vindija_Neanderthal == "1/0", "0/1", archaic_genotype_summ_Vindija_Neanderthal))
introgressed_variants_vindija_neanderthal_alt <- introgressed_variants_vindija_neanderthal %>% filter(is.na(VariantALT)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantCHROM, VariantPOS, archaic_genotype_summ_Vindija_Neanderthal)

introgressed_variants_chagyrskaya_neanderthal_ref <- introgressed_variants_chagyrskaya_neanderthal %>% filter(!is.na(VariantALT)) %>% 
	mutate(archaic_genotype_summ_Chagyrskaya_Neanderthal = ifelse(archaic_genotype_summ_Chagyrskaya_Neanderthal == "1/0", "0/1", archaic_genotype_summ_Chagyrskaya_Neanderthal))
introgressed_variants_chagyrskaya_neanderthal_alt <- introgressed_variants_chagyrskaya_neanderthal %>% filter(is.na(VariantALT)) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantCHROM, VariantPOS, archaic_genotype_summ_Chagyrskaya_Neanderthal)

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)

names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_altai_denisovan_ref))
introgressed_variants_altai_denisovan_ref <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_altai_denisovan_ref %>% mutate_at(names_temp, as.character))
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_altai_denisovan_alt))
introgressed_variants_altai_denisovan_alt <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_altai_denisovan_alt %>% mutate_at(names_temp, as.character))

introgressed_variants_altai_denisovan <- full_join(
	introgressed_variants_altai_denisovan_ref, 
	introgressed_variants_altai_denisovan_alt
)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_altai_denisovan))
introgressed_variants_altai_denisovan <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_altai_denisovan %>% mutate_at(names_temp, as.character))

names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_altai_neanderthal_ref))
introgressed_variants_altai_neanderthal_ref <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_altai_neanderthal_ref %>% mutate_at(names_temp, as.character))
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_altai_neanderthal_alt))
introgressed_variants_altai_neanderthal_alt <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_altai_neanderthal_alt %>% mutate_at(names_temp, as.character))

introgressed_variants_altai_neanderthal <- full_join(
	introgressed_variants_altai_neanderthal_ref, 
	introgressed_variants_altai_neanderthal_alt
)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_altai_neanderthal))
introgressed_variants_altai_neanderthal <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_altai_neanderthal %>% mutate_at(names_temp, as.character))

names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_vindija_neanderthal_ref))
introgressed_variants_vindija_neanderthal_ref <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_vindija_neanderthal_ref %>% mutate_at(names_temp, as.character))
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_vindija_neanderthal_alt))
introgressed_variants_vindija_neanderthal_alt <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_vindija_neanderthal_alt %>% mutate_at(names_temp, as.character))

introgressed_variants_vindija_neanderthal <- full_join(
	introgressed_variants_vindija_neanderthal_ref, 
	introgressed_variants_vindija_neanderthal_alt
)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_vindija_neanderthal))
introgressed_variants_vindija_neanderthal <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_vindija_neanderthal %>% mutate_at(names_temp, as.character))

names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_chagyrskaya_neanderthal_ref))
introgressed_variants_chagyrskaya_neanderthal_ref <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_chagyrskaya_neanderthal_ref %>% mutate_at(names_temp, as.character))
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_chagyrskaya_neanderthal_alt))
introgressed_variants_chagyrskaya_neanderthal_alt <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	inner_join(introgressed_variants_chagyrskaya_neanderthal_alt %>% mutate_at(names_temp, as.character))

introgressed_variants_chagyrskaya_neanderthal <- full_join(
	introgressed_variants_chagyrskaya_neanderthal_ref, 
	introgressed_variants_chagyrskaya_neanderthal_alt
)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_chagyrskaya_neanderthal))
introgressed_variants_chagyrskaya_neanderthal <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_chagyrskaya_neanderthal %>% mutate_at(names_temp, as.character))

# merge across genomes
introgressed_variants_archaic_genotypes <- introgressed_variants_altai_denisovan %>% 
	full_join(introgressed_variants_altai_neanderthal) %>% 
	full_join(introgressed_variants_vindija_neanderthal) %>% 
	full_join(introgressed_variants_chagyrskaya_neanderthal)

# save all variants
write_tsv(introgressed_variants_archaic_genotypes, gzfile("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_archaic_genotypes.txt.gz"))

# save summ variants
introgressed_variants_archaic_genotypes_summ <- introgressed_variants_archaic_genotypes %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_archaic_genotypes_summ, gzfile("../../results/2e-annotate_variants_archaic_genotypes/introgressed_variants_archaic_genotypes_summ.txt.gz"))
