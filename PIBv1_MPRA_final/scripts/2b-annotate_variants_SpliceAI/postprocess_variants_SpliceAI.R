#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)
library(genekitr)

# load VCF
introgressed_variants_SpliceAI <- as_tibble(read.vcfR("../../results/2b-annotate_variants_SpliceAI/introgressed_variants-SpliceAI_raw_hg19.vcf.gz")@fix) 

# add GRange and Variant cols
introgressed_variants_SpliceAI <- introgressed_variants_SpliceAI %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, INFO) 

# split INFO cols and add prefix
introgressed_variants_SpliceAI <- introgressed_variants_SpliceAI %>% 
	separate(col=INFO, into=paste("SpliceAI_raw_orig", c("ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"), sep="_"), remove=TRUE, sep="\\|", convert=TRUE) %>% 
	mutate(SpliceAI_raw_orig_ALLELE = gsub(".*=", "", SpliceAI_raw_orig_ALLELE)) 

# add temp SpliceAI MAX col
introgressed_variants_SpliceAI <- introgressed_variants_SpliceAI %>% 
	mutate(SpliceAI_raw_temp_DS_MAX = pmax(SpliceAI_raw_orig_DS_AG, SpliceAI_raw_orig_DS_AL, SpliceAI_raw_orig_DS_DG, SpliceAI_raw_orig_DS_DL))

# add temp gene id conversions
transId_map <- deframe(transId(sort(unique(introgressed_variants_SpliceAI$SpliceAI_raw_orig_SYMBOL)), transTo=c("symbol"), keepNA=TRUE, unique=TRUE, hgVersion="v38"))

introgressed_variants_SpliceAI <- introgressed_variants_SpliceAI %>% 
	rowwise() %>% mutate(SpliceAI_temp_Symbol = ifelse(SpliceAI_raw_orig_SYMBOL %in% names(transId_map), transId_map[[SpliceAI_raw_orig_SYMBOL]], NA)) %>% ungroup()

# add summary cols
introgressed_variants_SpliceAI <- introgressed_variants_SpliceAI %>% 
	mutate(SpliceAI_raw_summ_Gene = SpliceAI_raw_orig_SYMBOL) %>% 
	mutate(SpliceAI_raw_summ_Description = paste(SpliceAI_raw_orig_DS_AG, SpliceAI_raw_orig_DS_AL, SpliceAI_raw_orig_DS_DG, SpliceAI_raw_orig_DS_DL, SpliceAI_raw_orig_DP_AG, SpliceAI_raw_orig_DP_AL, SpliceAI_raw_orig_DP_DG, SpliceAI_raw_orig_DP_DL, sep="|")) %>% 
	mutate(SpliceAI_raw_summ_Score = SpliceAI_raw_temp_DS_MAX) %>% 
	mutate(SpliceAI_raw_summ_Summary = ifelse(SpliceAI_raw_temp_DS_MAX>=0.8, "SpliceAI>=0.8", ifelse(SpliceAI_raw_temp_DS_MAX>=0.5, "SpliceAI>=0.5", ifelse(SpliceAI_raw_temp_DS_MAX>=0.2, "SpliceAI>=0.2", NA)))) %>% 
	mutate(SpliceAI_raw_summ_Symbol = SpliceAI_temp_Symbol)

# merge to full variants
introgressed_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_temp <- intersect(names(introgressed_variants_tb), names(introgressed_variants_SpliceAI))
introgressed_variants_SpliceAI <- introgressed_variants_tb %>% mutate_at(names_temp, as.character)%>% 
	left_join(introgressed_variants_SpliceAI %>% mutate_at(names_temp, as.character))

# save all variants
write_tsv(introgressed_variants_SpliceAI, gzfile("../../results/2b-annotate_variants_SpliceAI/introgressed_variants_SpliceAI.txt.gz"))

# save summ variants
introgressed_variants_SpliceAI_summ <- introgressed_variants_SpliceAI %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(introgressed_variants_SpliceAI_summ, gzfile("../../results/2b-annotate_variants_SpliceAI/introgressed_variants_SpliceAI_summ.txt.gz"))
sort(table(introgressed_variants_SpliceAI_summ$SpliceAI_raw_summ_Summary))

# visual check
pdf("../../results/2b-annotate_variants_SpliceAI/introgressed_variants_SpliceAI_pie.pdf")
plot(sort(introgressed_variants_SpliceAI$SpliceAI_raw_summ_Score), ylab="SpliceAI max score")
dev.off()

pdf("../../results/2b-annotate_variants_SpliceAI/introgressed_variants_SpliceAI_summ_pie.pdf")
pie(rev(sort(table(introgressed_variants_SpliceAI_summ$SpliceAI_raw_summ_Summary))))
dev.off()
