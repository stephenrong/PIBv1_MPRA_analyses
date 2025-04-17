#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(wrapr)

# for bio data analysis
library(plyranges)
library(rtracklayer)
library(vcfR)
source("../shared_functions/seqinfo_fix_change.R")
source("../shared_functions/find_overlap_map.R")
source("../shared_functions/convert_to_vcf.R")

# list of pop names
pop_list <- c("Ata", "Baining", "Mamusi", "Melamela", "Lavongai")

# load variants
#   new version
adaptive_variants_tb <- as_tibble(fread("../../data/PIBv1_MPRA_pilot_extras/PIB92_pilot_MPRA_variants_Sprime_NEADENgt0.3_PIBfreqgt0.3_gSD_35mer_filtered_SNPs_noarcREF_nospanningDEL_wAlleles_match_VariantID_TractID.tsv.gz")) %>% 
	mutate(seqnames = as.character(seqnames), VariantCHROM = as.character(VariantCHROM))

# add population AFs
convert_vcf_to_tidy <- function(file_name, pop_name) {
	vcf_tb <- suppressWarnings(read.vcfR(file_name))
	vcf_tb <- bind_cols(
		vcf_tb@fix %>% as_tibble(), 
		vcf_tb %>% extract_info_tidy(info_fields = c("AC", "AF", "AN"), info_types = TRUE, info_sep = ";")
	) %>% 
		dplyr::select(-INFO, -Key)
	vcf_tb <- vcf_tb %>% 
		mutate(VariantCHROM = CHROM, VariantPOS = POS, VariantREF = REF, VariantALT = ALT) %>% 
		mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
		mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS, width=1, strand="*") %>% 
		dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, AF, AC, AN)
	names(vcf_tb)[11:13] <- paste(pop_name, names(vcf_tb)[11:13], sep="_")
	return(vcf_tb)
}

adaptive_variants_AFs <- 
	lapply(pop_list, function(pop_name) {
		convert_vcf_to_tidy(paste("../../data/PIBv1_MPRA_pilot_AFs/PIB92_SD_hs37d5_SNP_VQSR99.5_INDEL_VQSR99.0_dbSNP138_VQSRpass_nomissinggenos_autosomes_", 
			pop_name, "_filteredSamples_filteredSNPs_fillAFs_nogeno.vcf.gz", sep=""), pop_name)
	}) %>% 
	plyr::join_all(type="full") %>% as_tibble()

adaptive_variants_tb <- adaptive_variants_tb %>% 
	left_join(adaptive_variants_AFs %>% 
		dplyr::select(VariantID, matches("_AF|_AC|_AN")))

adaptive_variants_tb$PIB_AN <- rowSums(adaptive_variants_tb[,sapply(pop_list, function(x) {paste0(x, "_AN")})] %>% mutate_all(as.numeric))
adaptive_variants_tb$PIB_AC <- rowSums(adaptive_variants_tb[,sapply(pop_list, function(x) {paste0(x, "_AC")})] %>% mutate_all(as.numeric))
adaptive_variants_tb$PIB_AF <- adaptive_variants_tb$PIB_AC/adaptive_variants_tb$PIB_AN

# add other adaptive tracts
adaptive_variants_tracts <- as_tibble(fread("../../data/PIBv1_MPRA_pilot_extras/PIB92_pilot_MPRA_variants_Sprime_NEADENgt0.3_PIBfreqgt0.3_gSD_35mer_filtered_SNPs_noarcREF_nospanningDEL_wAlleles_filter_adaptive_TractID.tsv.gz")) %>% 
	mutate(seqnames = as.character(seqnames), VariantCHROM = as.character(VariantCHROM)) %>% 
	filter(VariantID %in% adaptive_variants_tb$VariantID) %>% 
	pivot_wider(names_from=Pop, values_from=TractID, names_glue="{Pop}_TractID", names_sort=TRUE)

adaptive_variants_tb <- adaptive_variants_tb %>% 
	left_join(adaptive_variants_tracts)

# replace adaptive tracts
# 	with other adaptive tracts
adaptive_tracts_tb <- NULL
for (pop in pop_list) {
	let(c(VAR=paste0(pop, "_TractID")), 
		adaptive_variants_tb_pop <- adaptive_variants_tb %>% 
			filter(!is.na(VAR)) %>% 
			mutate(Pop = pop)
	)
	let(c(VAR=paste0(pop, "_TractID")), 
		adaptive_variants_tb_pop_summ <- adaptive_variants_tb_pop %>% 
			group_by(seqnames, Pop, VAR, strand) %>% 
			summarise(
				start = min(start), 
				end = max(end)
			) %>% ungroup() %>% 
			dplyr::rename(TractID = VAR) %>% 
			dplyr::select(seqnames, start, end, strand, Pop, TractID)
	)
	if (is.null(adaptive_tracts_tb)) {
		adaptive_tracts_tb <- adaptive_variants_tb_pop_summ
	} else {
		adaptive_tracts_tb <- bind_rows(adaptive_tracts_tb, adaptive_variants_tb_pop_summ)
	}
}

# flag adaptive
adaptive_variants_tb <- adaptive_variants_tb %>% 
	mutate(Flag_variant_not_in_adaptive_tract = !(TractID %in% adaptive_tracts_tb$TractID))

# load oligos
mpra_oligos_tb <- as_tibble(fread("../../data/PIBv1_MPRA_pilot_oligos/PIB92_pilot_MPRA_final_oligos_inclControls_toOrder_SRb.txt.gz", header=F))
names(mpra_oligos_tb) <- c("Id", "Seq")
mpra_oligos_tb <- mpra_oligos_tb %>% filter(!grepl("_ctrl_", Id))

mpra_ccre_variants_tb <- mpra_oligos_tb %>% 
	separate(col=Id, into=c("VariantCHROM", "VariantPOS", "Allele", "Pop", "Sprime1", "Sprime2"), sep="_", remove=FALSE) %>% 
	mutate(TractID = paste(Pop, Sprime1, Sprime2, sep="_")) %>% dplyr::select(-Sprime1, -Sprime2) %>% 
	dplyr::select(Pop, TractID, VariantCHROM, VariantPOS, Allele, Pop, Id, Seq) %>% 
	mutate(VariantCHROM = gsub("chr", "", VariantCHROM), VariantPOS = as.numeric(VariantPOS))

mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	pivot_wider(names_from=Allele, values_from=c(Id, Seq)) %>% 
	mutate(VariantREF = substr(Seq_ref, 101, 101), VariantALT = substr(Seq_alt, 101, 101)) %>% 
	mutate(VariantID = paste(VariantCHROM, VariantPOS, VariantREF, VariantALT, sep="_")) %>% 
	dplyr::rename(IdRef=Id_ref, IdAlt=Id_alt, SeqRef=Seq_ref, SeqAlt=Seq_alt) %>% 
	dplyr::select(VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, Pop, TractID, IdRef, IdAlt, SeqRef, SeqAlt)

mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	mutate(seqnames=VariantCHROM, start=VariantPOS, end=VariantPOS+nchar(VariantREF)-1, width=nchar(VariantREF), strand="*") %>% 
	dplyr::select(seqnames, start, end, width, strand, VariantID, VariantCHROM, VariantPOS, VariantREF, VariantALT, Pop, TractID, IdRef, IdAlt, SeqRef, SeqAlt)

# join adaptive variant annotations
mpra_ccre_variants_tb <- adaptive_variants_tb %>% 
	right_join(mpra_ccre_variants_tb)

# flag liftOver
hg19ToHg38 <- import.chain("../../../Datasets/reference_genomes/liftOver_chains/hg19ToHg38.over.chain")

adaptive_variants_lift38_gr <- adaptive_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
mpra_ccre_variants_lift38_gr <- mpra_ccre_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()

adaptive_variants_tb <- adaptive_variants_tb %>% 
	mutate(Flag_variant_failed_liftOver_hg38 = !(VariantID %in% adaptive_variants_lift38_gr$VariantID))
mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	mutate(Flag_variant_failed_liftOver_hg38 = !(VariantID %in% mpra_ccre_variants_lift38_gr$VariantID))

# is in variant set
adaptive_variants_tb <- adaptive_variants_tb %>% 
	mutate(PIB_pilot_adaptive_variant = VariantID %in% adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_mpra_ccre_variant = (VariantID %in% mpra_ccre_variants_tb$VariantID))

mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	mutate(PIB_pilot_adaptive_variant = VariantID %in% adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_mpra_ccre_variant = (VariantID %in% mpra_ccre_variants_tb$VariantID))

# sort
adaptive_tracts_tb <- adaptive_tracts_tb %>% 
	GRanges() %>% sort() %>% as_tibble()
adaptive_variants_tb <- adaptive_variants_tb %>% 
	arrange(VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	GRanges() %>% sort() %>% as_tibble()
mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	arrange(VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	GRanges() %>% sort() %>% as_tibble()

# create GRanges
adaptive_tracts_gr <- adaptive_tracts_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()

adaptive_variants_gr <- adaptive_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
mpra_ccre_variants_gr <- mpra_ccre_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()

adaptive_variants_lift38_gr <- adaptive_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
mpra_ccre_variants_lift38_gr <- mpra_ccre_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()

# save as tables
write_tsv(adaptive_tracts_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_tracts.txt.gz"))
write_tsv(adaptive_variants_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz"))
write_tsv(mpra_ccre_variants_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_ccre_variants.txt.gz"))

# save as GRanges
saveRDS(adaptive_tracts_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_tracts.rds")
saveRDS(adaptive_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.rds")
saveRDS(mpra_ccre_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_ccre_variants.rds")

# save as BED files
rtracklayer::export(adaptive_tracts_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_tracts.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_ccre_variants.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_lift38_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38/adaptive_variants_lift38.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_lift38_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38/mpra_ccre_variants_lift38.bed.gz", format="BED")

# save as BED with chr
rtracklayer::export(adaptive_tracts_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_ranges_chr/adaptive_tracts_chr.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_chr/adaptive_variants_chr.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_chr/mpra_ccre_variants_chr.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_lift38_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_lift38_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/mpra_ccre_variants_lift38_chr.bed.gz", format="BED")

# save as VCF files
convert_gr_to_vcf(adaptive_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.vcf")
convert_gr_to_vcf(mpra_ccre_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_ccre_variants.vcf")
convert_gr_to_vcf(adaptive_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/mpra_ccre_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38/adaptive_variants_lift38.vcf")
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38/mpra_ccre_variants_lift38.vcf")
convert_gr_to_vcf(adaptive_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38/adaptive_variants_lift38_clean.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38/mpra_ccre_variants_lift38_clean.vcf", clean=TRUE)

# save as VCF with chr
convert_gr_to_vcf(adaptive_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_chr/adaptive_variants_chr.vcf")
convert_gr_to_vcf(mpra_ccre_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_chr/mpra_ccre_variants_chr.vcf")
convert_gr_to_vcf(adaptive_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_chr/adaptive_variants_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_chr/mpra_ccre_variants_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.vcf")
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/mpra_ccre_variants_lift38_chr.vcf")
convert_gr_to_vcf(adaptive_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/adaptive_variants_lift38_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_pilot/alternate_variants_lift38_chr/mpra_ccre_variants_lift38_clean_chr.vcf", clean=TRUE)

# save flagged
adaptive_variants_flagged_adaptive_tb <- filter(adaptive_variants_tb, Flag_variant_not_in_adaptive_tract)
adaptive_variants_flagged_liftOver_tb <- filter(adaptive_variants_tb, Flag_variant_failed_liftOver_hg38)
mpra_ccre_variants_flagged_adaptive_tb <- filter(mpra_ccre_variants_tb, Flag_variant_not_in_adaptive_tract)
mpra_ccre_variants_flagged_liftOver_tb <- filter(mpra_ccre_variants_tb, Flag_variant_failed_liftOver_hg38)
write_tsv(adaptive_variants_flagged_adaptive_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/flagged_variants/adaptive_variants_flagged_adaptive.txt.gz"))
write_tsv(adaptive_variants_flagged_liftOver_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/flagged_variants/adaptive_variants_flagged_liftOver.txt.gz"))
write_tsv(mpra_ccre_variants_flagged_adaptive_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/flagged_variants/mpra_ccre_variants_flagged_adaptive.txt.gz"))
write_tsv(mpra_ccre_variants_flagged_liftOver_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_pilot/flagged_variants/mpra_ccre_variants_flagged_liftOver.txt.gz"))
