#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)

# for bio data analysis
library(plyranges)
library(vcfR)
library(genekitr)

# ensembl vep order
adaptive_variants_Ensembl_VEP_order <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", "splice_region_variant", "splice_donor_5th_base_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant", "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation", "regulatory_region_variant", "feature_truncation", "intergenic_variant")

# ensembl vep names
adaptive_variants_Ensembl_VEP_names <- c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "Extra", "IMPACT", "SYMBOL", "SYMBOL_SOURCE", "STRAND", "BIOTYPE")

# load VCF 
adaptive_variants_Ensembl_VEP <- as_tibble(read.vcfR("../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants-Ensembl_VEP_nopick.vcf.gz")@fix) 

# add GRange and Variant cols
adaptive_variants_Ensembl_VEP <- adaptive_variants_Ensembl_VEP %>% 
	mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, INFO)

# split INFO column and add prefix
adaptive_variants_Ensembl_VEP <- adaptive_variants_Ensembl_VEP %>% 
	mutate(INFO = gsub("CSQ=", "", INFO)) %>% 
	separate_rows(INFO, sep=",", convert=TRUE) %>% 
	separate(col=INFO, into=paste("Ensembl_VEP_orig", adaptive_variants_Ensembl_VEP_names, sep="_"), remove=TRUE, sep="\\|") 

# add temp Consequence cols
adaptive_variants_Ensembl_VEP <- adaptive_variants_Ensembl_VEP %>% 
	mutate(Ensembl_VEP_temp_Consequence_no_modifier = gsub("&.*", "", Ensembl_VEP_orig_Consequence)) %>% 
	mutate(Ensembl_VEP_temp_Consequence_impact_only = ifelse((Ensembl_VEP_orig_IMPACT %in% c("HIGH", "MODERATE") | Ensembl_VEP_orig_Consequence %in% c("splice_region_variant", "5_prime_UTR_variant", "3_prime_UTR_variant")), Ensembl_VEP_temp_Consequence_no_modifier, NA))

# add temp gene id conversions
transId_map <- deframe(transId(sort(unique(adaptive_variants_Ensembl_VEP$Ensembl_VEP_orig_Gene)), transTo=c("symbol"), keepNA=TRUE, unique=TRUE, hgVersion="v38"))
adaptive_variants_Ensembl_VEP <- adaptive_variants_Ensembl_VEP %>% 
	rowwise() %>% mutate(Ensembl_VEP_temp_Symbol = ifelse(Ensembl_VEP_orig_Gene != "", transId_map[[Ensembl_VEP_orig_Gene]], NA)) %>% ungroup()

# save names
adaptive_variants_Ensembl_VEP_names <- names(adaptive_variants_Ensembl_VEP)[11:length(adaptive_variants_Ensembl_VEP)]

# collapse multiple entries
adaptive_variants_Ensembl_VEP <- adaptive_variants_Ensembl_VEP %>% 
	group_by(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	arrange(match(Ensembl_VEP_temp_Consequence_no_modifier, adaptive_variants_Ensembl_VEP_order)) %>% 
	dplyr::summarise_at(adaptive_variants_Ensembl_VEP_names, function(x) {paste(x, collapse=",")}) %>% 
	ungroup()

# add summary cols
adaptive_variants_Ensembl_VEP <- adaptive_variants_Ensembl_VEP %>% 
	mutate(Ensembl_VEP_summ_Gene = gsub(",.*", "", Ensembl_VEP_orig_Gene)) %>% 
	mutate(Ensembl_VEP_summ_Gene = ifelse(Ensembl_VEP_summ_Gene != "", Ensembl_VEP_summ_Gene, NA)) %>% 
	mutate(Ensembl_VEP_summ_Category = gsub(",.*", "", gsub(",NA|^NA", "", Ensembl_VEP_temp_Consequence_no_modifier))) %>% 
	mutate(Ensembl_VEP_summ_Summary = gsub(",.*", "", gsub(",NA|^NA", "", Ensembl_VEP_temp_Consequence_impact_only))) %>% 
	mutate(Ensembl_VEP_summ_Summary = ifelse(Ensembl_VEP_summ_Summary != "", Ensembl_VEP_summ_Summary, NA)) %>% 
	mutate(Ensembl_VEP_summ_Symbol = gsub(",.*", "", Ensembl_VEP_temp_Symbol)) %>% 
	mutate(Ensembl_VEP_summ_Symbol = ifelse(Ensembl_VEP_summ_Symbol != "NA", Ensembl_VEP_summ_Symbol, NA))

# merge to full variants
adaptive_variants_tb <- as_tibble(fread("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz")) %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT)
names_inter <- intersect(names(adaptive_variants_tb), names(adaptive_variants_Ensembl_VEP))
adaptive_variants_Ensembl_VEP <- adaptive_variants_tb %>% mutate_at(names_inter, as.character) %>% 
	left_join(adaptive_variants_Ensembl_VEP %>% mutate_at(names_inter, as.character))

# save all variants
write_tsv(adaptive_variants_Ensembl_VEP, gzfile("../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants_Ensembl_VEP.txt.gz"))

# save summ variants
adaptive_variants_Ensembl_VEP_summ <- adaptive_variants_Ensembl_VEP %>% 
	dplyr::select(seqnames, start, end, width, strand, starts_with("Variant"), contains("_summ_"))
write_tsv(adaptive_variants_Ensembl_VEP_summ, gzfile("../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants_Ensembl_VEP_summ.txt.gz"))
sort(table(adaptive_variants_Ensembl_VEP_summ$Ensembl_VEP_summ_Summary))

# visual checks
pdf("../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants_Ensembl_VEP_pie.pdf")
pie(rev(sort(table(adaptive_variants_Ensembl_VEP$Ensembl_VEP_summ_Category))))
dev.off()

pdf("../../results/2a-annotate_variants_Ensembl_VEP/adaptive_variants_Ensembl_VEP_summ_pie.pdf")
pie(rev(sort(table(adaptive_variants_Ensembl_VEP_summ$Ensembl_VEP_summ_Summary))))
dev.off()
