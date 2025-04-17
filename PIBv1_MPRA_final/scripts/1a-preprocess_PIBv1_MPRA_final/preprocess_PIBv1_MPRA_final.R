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

map_pop_to_pilot_pop <- c(
	"Ata"="Ata", 
	"Baining-Kagat"="Baining", 
	"Baining-Mali"="Baining", 
	"Lavongai-Mussau"="Lavongai", 
	"Mamusi"="Mamusi", 
	"Melamela"="Melamela"
)

map_saf_to_col_name <- c(
	"Ata"="Ata", 
	"BainingKagat"="BainingKagat", 
	"BainingMali"="BainingMali", 
	"Bellona"="BellonaRennell", 
	"Goroka"="GorokaSepik", 
	"Kove"="Kove", 
	"LavongaiMussau"="LavongaiMussau", 
	"Malaita"="Malaita", 
	"Mamusi"="Mamusi", 
	"Melamela"="Melamela", 
	"NailikNotsiTigak"="NailikNotsiTigak", 
	"Nakanai"="NakanaiMangseng", 
	"Nasioi"="Nasioi", 
	"Rennell"=NA, 
	"SantaCruz"="SantaCruz", 
	"Saposa"="Saposa", 
	"Sepik"=NA, 
	"Tikopia"="Tikopia", 
	"VellaLavella"="VellaLavella"
)

misc_saf_list <- c(
	"AFR", "AMR", "CSA", "EAS", "EUR", "ISEA", "MDE", "OCN"
)

map_pop_to_col_name_short <- map_pop_to_col_name[which(!is.na(map_pop_to_col_name))]
map_saf_to_col_name_short <- map_saf_to_col_name[which(!is.na(map_saf_to_col_name))]

# load introgressed tracts
introgressed_tracts_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Sprime_tracts_OCNnoVanuatu_MRfilter_wGeneList.bed")) %>% 
	# convert from 0-indexed to 1-indexed
	mutate(Start = Start + 1)

# add classifications
introgressed_tracts_nean <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Neandertal_Sprime_tracts_OCNnoVanuatu_MRfilter_wGeneList.bed")) %>% 
	mutate(Start = Start + 1)
introgressed_tracts_deni <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Denisovan_Sprime_tracts_OCNnoVanuatu_MRfilter_wGeneList.bed")) %>% 
	mutate(Start = Start + 1)
introgressed_tracts_ambig <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Ambiguous_Sprime_tracts_OCNnoVanuatu_MRfilter_wGeneList.bed")) %>% 
	mutate(Start = Start + 1)

introgressed_tracts_tb <- introgressed_tracts_tb %>% 
	mutate(NeandertalClass = (TractID %in% introgressed_tracts_nean$TractID)) %>% 
	mutate(DenisovanClass = (TractID %in% introgressed_tracts_deni$TractID)) %>% 
	mutate(AmbiguousClass = (TractID %in% introgressed_tracts_ambig$TractID))

# drop extra pops
introgressed_tracts_tb <- introgressed_tracts_tb %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# 	fix columns
names(introgressed_tracts_tb)[4:18] <- paste0("Tracts_", names(introgressed_tracts_tb)[4:18])
introgressed_tracts_tb <- introgressed_tracts_tb %>% 
	mutate(
		seqnames = `#Chromosome`, 
		start = Start,
		end = End, 
		width = End - Start + 1, 
		strand = "*"
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Tracts_")
	)

# load introgressed tracts out
introgressed_tracts_noOCN_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed")) %>% 
	# convert from 0-indexed to 1-indexed
	mutate(Start = Start + 1)

# add classifications
introgressed_tracts_noOCN_nean <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Neandertal_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed")) %>% 
	mutate(Start = Start + 1)
introgressed_tracts_noOCN_deni <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Denisovan_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed")) %>% 
	mutate(Start = Start + 1)
introgressed_tracts_noOCN_ambig <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Ambiguous_Sprime_tracts_noVanuatu_MRfilter_wGeneList.bed")) %>% 
	mutate(Start = Start + 1)

introgressed_tracts_noOCN_tb <- introgressed_tracts_noOCN_tb %>% 
	mutate(NeandertalClass = (TractID %in% introgressed_tracts_noOCN_nean$TractID)) %>% 
	mutate(DenisovanClass = (TractID %in% introgressed_tracts_noOCN_deni$TractID)) %>% 
	mutate(AmbiguousClass = (TractID %in% introgressed_tracts_noOCN_ambig$TractID))

# 	fix columns
names(introgressed_tracts_noOCN_tb)[4:18] <- paste0("Tracts_", names(introgressed_tracts_noOCN_tb)[4:18])
introgressed_tracts_noOCN_tb <- introgressed_tracts_noOCN_tb %>% 
	mutate(
		seqnames = `#Chromosome`, 
		start = Start,
		end = End, 
		width = End - Start + 1, 
		strand = "*"
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Tracts_")
	)

# 	remove tracts
introgressed_tracts_noOCN_tb <- introgressed_tracts_noOCN_tb %>% 
	filter(!(Tracts_SprimePopulation %in% names(map_pop_to_col_name)))

# load adaptive tracts
adaptive_tracts_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Sprime_putative_adaptive_introgressed_tracts_OCNnoVanuatu_MRfilterBeforeTFfilter_wGeneList.bed")) %>% 
	# convert from 0-indexed to 1-indexed
	mutate(Start = Start + 1)

# add classifications
adaptive_tracts_tb <- adaptive_tracts_tb %>% 
	mutate(NeandertalClass = (TractID %in% introgressed_tracts_nean$TractID)) %>% 
	mutate(DenisovanClass = (TractID %in% introgressed_tracts_deni$TractID)) %>% 
	mutate(AmbiguousClass = (TractID %in% introgressed_tracts_ambig$TractID))

# drop extra pops
adaptive_tracts_tb <- adaptive_tracts_tb %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# 	fix columns
names(adaptive_tracts_tb)[4:18] <- paste0("Tracts_", names(adaptive_tracts_tb)[4:18])
adaptive_tracts_tb <- adaptive_tracts_tb %>% 
	mutate(
		seqnames = `#Chromosome`, 
		start = Start,
		end = End,
		width = End - Start + 1, 
		strand = "*"
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Tracts_")
	)

# load adaptive corehaps
adaptive_corehaps_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Sprime_putative_adaptive_introgressed_corehaps_OCNnoVanuatu_wGeneList.bed")) %>% 
	# convert from 0-indexed to 1-indexed
	mutate(Start = Start + 1)

# drop extra pops
adaptive_corehaps_tb <- adaptive_corehaps_tb %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# 	fix columns
names(adaptive_corehaps_tb)[4:10] <- paste0("Corehaps_", names(adaptive_corehaps_tb)[4:10])
adaptive_corehaps_tb <- adaptive_corehaps_tb %>% 
	mutate(
		seqnames = `#Chromosome`, 
		start = Start,
		end = End, 
		width = End - Start + 1, 
		strand = "*"
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Corehaps_")
	)

# load introgressed corehaps
introgressed_corehaps_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_Sprime_corehaps_OCNnoVanuatu_wGeneList.bed")) %>% 
	# convert from 0-indexed to 1-indexed
	mutate(Start = Start + 1)

# drop extra pops
introgressed_corehaps_tb <- introgressed_corehaps_tb %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# 	fix columns
names(introgressed_corehaps_tb)[4:10] <- paste0("Corehaps_", names(introgressed_corehaps_tb)[4:10])
introgressed_corehaps_tb <- introgressed_corehaps_tb %>% 
	mutate(
		seqnames = `#Chromosome`, 
		start = Start,
		end = End, 
		width = End - Start + 1, 
		strand = "*"
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Corehaps_")
	)

# load introgressed variants SAF
introgressed_variants_SAF_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_OCNnoVanuatu_SprimeAlleles_alltracts_wSprimeAFs.tsv.gz"))

introgressed_variants_SAF_nean <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_OCNnoVanuatu_SprimeAlleles_Neandertaltracts_wSprimeAFs.tsv.gz"))
introgressed_variants_SAF_deni <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_OCNnoVanuatu_SprimeAlleles_Denisovantracts_wSprimeAFs.tsv.gz"))
introgressed_variants_SAF_ambig <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_OCNnoVanuatu_SprimeAlleles_Ambiguoustracts_wSprimeAFs.tsv.gz"))

# drop extra pops
introgressed_variants_SAF_tb <- introgressed_variants_SAF_tb %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# normalize ALTs
introgressed_variants_SAF_tb <- introgressed_variants_SAF_tb %>% 
	separate_rows(ALT, sep=",")

# fix columns
names(introgressed_variants_SAF_tb)[5:15] <- paste0("Alleles_", names(introgressed_variants_SAF_tb)[5:15])
introgressed_variants_SAF_tb <- introgressed_variants_SAF_tb %>% 
	# add GRange columns
	mutate(
		seqnames = CHROM,
		start = POS,
		end = POS,
		width = 1,
		strand = "*"
	) %>% 
	# standardize variants
	mutate(
		VariantID = paste(CHROM, POS, REF, ALT, sep="_"),
		VariantCHROM = CHROM,
		VariantPOS = POS,
		VariantREF = REF,
		VariantALT = ALT
	)

# drop extra pops
introgressed_variants_SAF_nean <- introgressed_variants_SAF_nean %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# normalize ALTs
introgressed_variants_SAF_nean <- introgressed_variants_SAF_nean %>% 
	separate_rows(ALT, sep=",")

# fix columns
names(introgressed_variants_SAF_nean)[5:15] <- paste0("Alleles_", names(introgressed_variants_SAF_nean)[5:15])
introgressed_variants_SAF_nean <- introgressed_variants_SAF_nean %>% 
	# add GRange columns
	mutate(
		seqnames = CHROM,
		start = POS,
		end = POS,
		width = 1,
		strand = "*"
	) %>% 
	# standardize variants
	mutate(
		VariantID = paste(CHROM, POS, REF, ALT, sep="_"),
		VariantCHROM = CHROM,
		VariantPOS = POS,
		VariantREF = REF,
		VariantALT = ALT
	)

# drop extra pops
introgressed_variants_SAF_deni <- introgressed_variants_SAF_deni %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# normalize ALTs
introgressed_variants_SAF_deni <- introgressed_variants_SAF_deni %>% 
	separate_rows(ALT, sep=",")

# fix columns
names(introgressed_variants_SAF_deni)[5:15] <- paste0("Alleles_", names(introgressed_variants_SAF_deni)[5:15])
introgressed_variants_SAF_deni <- introgressed_variants_SAF_deni %>% 
	# add GRange columns
	mutate(
		seqnames = CHROM,
		start = POS,
		end = POS,
		width = 1,
		strand = "*"
	) %>% 
	# standardize variants
	mutate(
		VariantID = paste(CHROM, POS, REF, ALT, sep="_"),
		VariantCHROM = CHROM,
		VariantPOS = POS,
		VariantREF = REF,
		VariantALT = ALT
	)

# drop extra pops
introgressed_variants_SAF_ambig <- introgressed_variants_SAF_ambig %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# normalize ALTs
introgressed_variants_SAF_ambig <- introgressed_variants_SAF_ambig %>% 
	separate_rows(ALT, sep=",")

# fix columns
names(introgressed_variants_SAF_ambig)[5:15] <- paste0("Alleles_", names(introgressed_variants_SAF_ambig)[5:15])
introgressed_variants_SAF_ambig <- introgressed_variants_SAF_ambig %>% 
	# add GRange columns
	mutate(
		seqnames = CHROM,
		start = POS,
		end = POS,
		width = 1,
		strand = "*"
	) %>% 
	# standardize variants
	mutate(
		VariantID = paste(CHROM, POS, REF, ALT, sep="_"),
		VariantCHROM = CHROM,
		VariantPOS = POS,
		VariantREF = REF,
		VariantALT = ALT
	)

# add classifications
introgressed_variants_SAF_tb <- introgressed_variants_SAF_tb %>% 
	mutate(Alleles_NeandertalClass = (VariantID %in% introgressed_variants_SAF_nean$VariantID)) %>% 
	mutate(Alleles_DenisovanClass = (VariantID %in% introgressed_variants_SAF_deni$VariantID)) %>% 
	mutate(Alleles_AmbiguousClass = (VariantID %in% introgressed_variants_SAF_ambig$VariantID))

# pivot tracts by pop
introgressed_variants_SAF_tb_tractid <- introgressed_variants_SAF_tb %>% 
	dplyr::select(-Alleles_CorehapID) %>% 
	pivot_wider(
		names_from = Alleles_SprimePopulation, 
		names_glue = "Introgr_TractID_{Alleles_SprimePopulation}",
		values_from = Alleles_TractID
	) 

introgressed_variants_SAF_tb_corehapid <- introgressed_variants_SAF_tb %>% 
	dplyr::select(-Alleles_TractID) %>% 
	pivot_wider(
		names_from = Alleles_SprimePopulation, 
		names_glue = "Introgr_CorehapID_{Alleles_SprimePopulation}",
		values_from = Alleles_CorehapID
	) 

introgressed_variants_SAF_tb <- full_join(
	introgressed_variants_SAF_tb_tractid, 
	introgressed_variants_SAF_tb_corehapid
)

# clean up columns
introgressed_variants_SAF_tb <- introgressed_variants_SAF_tb %>% 
	dplyr::select(
		-CHROM, -POS, -REF, -ALT
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Variant"), 
		starts_with("Alleles_"), 
		starts_with("Introgr_TractID_"),
		starts_with("Introgr_CorehapID_"),
		starts_with("SAF_")
	)

# load introgressed variants AF
introgressed_variants_AF_tb <- as_tibble(fread("../../data/PIBv1_MPRA_final_extras/PIBv1_OCNnoVanuatu_SprimeAlleles_alltracts_wAFs.tsv.gz"))

# drop extra pops
introgressed_variants_AF_tb <- introgressed_variants_AF_tb %>% 
	filter(!is.na(map_pop_to_col_name[SprimePopulation]))

# normalize ALTs
introgressed_variants_AF_tb <- introgressed_variants_AF_tb %>% 
	separate_rows(ALT, sep=",")

# fix columns
names(introgressed_variants_AF_tb)[5:15] <- paste0("Alleles_", names(introgressed_variants_AF_tb)[5:15])
introgressed_variants_AF_tb <- introgressed_variants_AF_tb %>% 
	# add GRange columns
	mutate(
		seqnames = CHROM,
		start = POS,
		end = POS,
		width = 1,
		strand = "*"
	) %>% 
	# standardize variants
	mutate(
		VariantID = paste(CHROM, POS, REF, ALT, sep="_"),
		VariantCHROM = CHROM,
		VariantPOS = POS,
		VariantREF = REF,
		VariantALT = ALT
	)

# pivot tracts by pop
introgressed_variants_AF_tb_tractid <- introgressed_variants_AF_tb %>% 
	dplyr::select(-Alleles_CorehapID) %>% 
	pivot_wider(
		names_from = Alleles_SprimePopulation, 
		names_glue = "Introgr_TractID_{Alleles_SprimePopulation}",
		values_from = Alleles_TractID
	)

introgressed_variants_AF_tb_corehapid <- introgressed_variants_AF_tb %>% 
	dplyr::select(-Alleles_TractID) %>% 
	pivot_wider(
		names_from = Alleles_SprimePopulation, 
		names_glue = "Introgr_CorehapID_{Alleles_SprimePopulation}",
		values_from = Alleles_CorehapID
	) 

introgressed_variants_AF_tb <- full_join(
	introgressed_variants_AF_tb_tractid, 
	introgressed_variants_AF_tb_corehapid
)

# clean up columns
introgressed_variants_AF_tb <- introgressed_variants_AF_tb %>% 
	dplyr::select(
		-CHROM, -POS, -REF, -ALT
	) %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Variant"), 
		starts_with("Alleles_"), 
		starts_with("Introgr_TractID_"),
		starts_with("Introgr_CorehapID_"),
		starts_with("AF_")
	)

# join adaptive variants SAF and AF
introgressed_variants_tb <- full_join(introgressed_variants_SAF_tb, introgressed_variants_AF_tb)

# SprimeAllele should be either REF or ALT
# this is a byproduct of the normalization step
introgressed_variants_tb <- introgressed_variants_tb %>% 
	filter((Alleles_SprimeAllele == VariantREF) | (Alleles_SprimeAllele == VariantALT))

# SprimeAllele should not be both REF and ALT
# this may be biologically meaningful if for example there's ILS
introgressed_variants_tb <- introgressed_variants_tb %>% 
	group_by(VariantID) %>% filter(n() == 1) %>% ungroup()

# clean column names
names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("Introgr_TractID_", names(map_pop_to_col_name)))] <- 
	sapply(names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("Introgr_TractID_", names(map_pop_to_col_name)))], 
		function(x) {
			for (key in names(map_pop_to_col_name)) {
				x <- gsub(paste0("_", key, "$"), paste0("_", map_pop_to_col_name[key]), x)
			}
		return(x)
		}, USE.NAMES = FALSE)

names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("Introgr_CorehapID_", names(map_pop_to_col_name)))] <- 
	sapply(names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("Introgr_CorehapID_", names(map_pop_to_col_name)))], 
		function(x) {
			for (key in names(map_pop_to_col_name)) {
				x <- gsub(paste0("_", key, "$"), paste0("_", map_pop_to_col_name[key]), x)
			}
		return(x)
		}, USE.NAMES = FALSE)

# added to deal with SAF/AF alignment issue
introgressed_variants_tb <- introgressed_variants_tb %>% 
	dplyr::select(-paste0("SAF_", names(map_saf_to_col_name)[which(is.na(map_saf_to_col_name))]))
names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("SAF_", names(map_saf_to_col_name)))] <- 
	sapply(names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("SAF_", names(map_saf_to_col_name)))], 
		function(x) {
			for (key in names(map_saf_to_col_name)) {
				x <- gsub(paste0("_", key, "$"), paste0("_", map_saf_to_col_name[key]), x)
			}
		return(x)
		}, USE.NAMES = FALSE)

introgressed_variants_tb <- introgressed_variants_tb %>% 
	dplyr::select(-paste0("AF_", names(map_saf_to_col_name)[which(is.na(map_saf_to_col_name))]))
names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("AF_", names(map_saf_to_col_name)))] <- 
	sapply(names(introgressed_variants_tb)[which(names(introgressed_variants_tb) %in% paste0("AF_", names(map_saf_to_col_name)))], 
		function(x) {
			for (key in names(map_saf_to_col_name)) {
				x <- gsub(paste0("_", key, "$"), paste0("_", map_saf_to_col_name[key]), x)
			}
		return(x)
		}, USE.NAMES = FALSE)

# get adaptive tracts
for (i in map_pop_to_col_name_short) {
	let(c(VAR1 = paste0("Introgr_TractID_", i), VAR2 = paste0("TractID_", i)), 
		introgressed_variants_tb <- introgressed_variants_tb %>% 
			mutate(VAR2 = ifelse(VAR1 %in% adaptive_tracts_tb$Tracts_TractID, VAR1, NA))
	)
	let(c(VAR1 = paste0("Introgr_CorehapID_", i), VAR2 = paste0("CorehapID_", i)), 
		introgressed_variants_tb <- introgressed_variants_tb %>% 
			mutate(VAR2 = ifelse(VAR1 %in% adaptive_corehaps_tb$Corehaps_TractID, VAR1, NA))
	)
}

# add overlap counts
TractID_pop_list <- paste0("TractID_", map_pop_to_col_name_short)
CorehapID_pop_list <- paste0("CorehapID_", map_pop_to_col_name_short)
introgressed_variants_tb <- introgressed_variants_tb %>% 
	mutate(
		Tracts_count = rowSums(!is.na(across(all_of(TractID_pop_list)))),
		Corehaps_count = rowSums(!is.na(across(all_of(CorehapID_pop_list)))),
	)

# reorder all columns
introgressed_variants_tb <- introgressed_variants_tb %>% 
	dplyr::select(
		seqnames, start, end, width, strand, 
		starts_with("Variant"), 
		starts_with("Alleles_"), 
		paste0("Introgr_TractID_", map_pop_to_col_name_short), 
		paste0("Introgr_CorehapID_", map_pop_to_col_name_short), 
		paste0("TractID_", map_pop_to_col_name_short), 
		paste0("CorehapID_", map_pop_to_col_name_short), 
		Tracts_count, Corehaps_count,
		paste0("SAF_", map_saf_to_col_name_short), 
		paste0("AF_", map_saf_to_col_name_short),
		paste0("SAF_", misc_saf_list), 
		paste0("AF_", misc_saf_list)
	)

# filter to only SNPs
introgressed_variants_tb <- introgressed_variants_tb %>% 
	filter(
		VariantREF %in% c("A", "C", "T", "G"),
		VariantALT %in% c("A", "C", "T", "G")
	)

# get adaptive variants
adaptive_variants_tb <- introgressed_variants_tb %>% 
	filter(Tracts_count > 0)
corehaps_variants_tb <- introgressed_variants_tb %>% 
	filter(Corehaps_count > 0)

# overlap ENCODE cCREs
ENCODE_cCREs_tb <- as_tibble(fread("../../../Datasets/gene_regulation_element_catalogs/ENCODE_SCREEN_V4_cCREs/data_cleanup/sample_agnostic/ENCODE_SCREEN_V4_cCREs_lift37/GRCh38-cCREs.V4.lift37.nochr.sort.bed.gz"))
ENCODE_cCREs_tb <- ENCODE_cCREs_tb[,c(1,2,3,6)]
ENCODE_cCREs_tb <- unique(ENCODE_cCREs_tb)
names(ENCODE_cCREs_tb) <- c("seqnames", "start", "end", "ENCODE_cCREs_orig_class")
ENCODE_cCREs_tb <- ENCODE_cCREs_tb %>% 
	mutate(start = start + 1) %>%  # translate 0-indexed to 1-indexed
	mutate(ENCODE_cCREs_orig_element = paste(paste(seqnames, start, end, sep="_"), ENCODE_cCREs_orig_class, sep="-"))
mpra_ccre_variants_tb <- subsetByOverlaps(
	adaptive_variants_tb %>% GRanges(), 
	ENCODE_cCREs_tb %>% GRanges()
) %>% as_tibble()

# flag liftOver
hg19ToHg38 <- import.chain("../../../Datasets/reference_genomes/liftOver_chains/hg19ToHg38.over.chain")

introgressed_variants_lift38_gr <- introgressed_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
adaptive_variants_lift38_gr <- adaptive_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
corehaps_variants_lift38_gr <- corehaps_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
mpra_ccre_variants_lift38_gr <- mpra_ccre_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()

introgressed_variants_tb <- introgressed_variants_tb %>% 
	mutate(Flag_variant_failed_liftOver_hg38 = !(VariantID %in% introgressed_variants_lift38_gr$VariantID))
adaptive_variants_tb <- adaptive_variants_tb %>% 
	mutate(Flag_variant_failed_liftOver_hg38 = !(VariantID %in% adaptive_variants_lift38_gr$VariantID))
corehaps_variants_tb <- corehaps_variants_tb %>% 
	mutate(Flag_variant_failed_liftOver_hg38 = !(VariantID %in% corehaps_variants_lift38_gr$VariantID))
mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	mutate(Flag_variant_failed_liftOver_hg38 = !(VariantID %in% mpra_ccre_variants_lift38_gr$VariantID))

# is in variant set
pilot_adaptive_variants_tb <- as_tibble(fread("../../../PIBv1_MPRA_pilot/results/1a-preprocess_PIBv1_MPRA_pilot/adaptive_variants.txt.gz"))
pilot_mpra_ccre_variants_tb <- as_tibble(fread("../../../PIBv1_MPRA_pilot/results/1a-preprocess_PIBv1_MPRA_pilot/mpra_ccre_variants.txt.gz"))

map_pop_to_col_name[names(map_pop_to_pilot_pop)]

introgressed_variants_tb <- introgressed_variants_tb %>% 
	mutate(PIB_final_introgressed_variant = VariantID %in% introgressed_variants_tb$VariantID) %>% 
	mutate(PIB_final_adaptive_variant = VariantID %in% adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_final_corehaps_variant = VariantID %in% corehaps_variants_tb$VariantID) %>% 
	mutate(PIB_final_mpra_ccre_variant = VariantID %in% mpra_ccre_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_adaptive_variant = VariantID %in% pilot_adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_mpra_ccre_variant = (VariantID %in% pilot_mpra_ccre_variants_tb$VariantID)) %>% 
	mutate(PIB_final_in_pilot_adaptive_variant = 
		(rowSums(!is.na(across(all_of(paste0("TractID_", map_pop_to_col_name[names(map_pop_to_pilot_pop)]))))) > 0)) %>% 
	mutate(PIB_final_in_pilot_mpra_ccre_variant = (PIB_final_in_pilot_adaptive_variant & PIB_final_mpra_ccre_variant))

adaptive_variants_tb <- adaptive_variants_tb %>% 
	mutate(PIB_final_introgressed_variant = VariantID %in% introgressed_variants_tb$VariantID) %>% 
	mutate(PIB_final_adaptive_variant = VariantID %in% adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_final_corehaps_variant = VariantID %in% corehaps_variants_tb$VariantID) %>% 
	mutate(PIB_final_mpra_ccre_variant = VariantID %in% mpra_ccre_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_adaptive_variant = VariantID %in% pilot_adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_mpra_ccre_variant = (VariantID %in% pilot_mpra_ccre_variants_tb$VariantID)) %>% 
	mutate(PIB_final_in_pilot_adaptive_variant = 
		(rowSums(!is.na(across(all_of(paste0("TractID_", map_pop_to_col_name[names(map_pop_to_pilot_pop)]))))) > 0)) %>% 
	mutate(PIB_final_in_pilot_mpra_ccre_variant = (PIB_final_in_pilot_adaptive_variant & PIB_final_mpra_ccre_variant))

corehaps_variants_tb <- corehaps_variants_tb %>% 
	mutate(PIB_final_introgressed_variant = VariantID %in% introgressed_variants_tb$VariantID) %>% 
	mutate(PIB_final_adaptive_variant = VariantID %in% adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_final_corehaps_variant = VariantID %in% corehaps_variants_tb$VariantID) %>% 
	mutate(PIB_final_mpra_ccre_variant = VariantID %in% mpra_ccre_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_adaptive_variant = VariantID %in% pilot_adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_mpra_ccre_variant = (VariantID %in% pilot_mpra_ccre_variants_tb$VariantID)) %>% 
	mutate(PIB_final_in_pilot_adaptive_variant = 
		(rowSums(!is.na(across(all_of(paste0("TractID_", map_pop_to_col_name[names(map_pop_to_pilot_pop)]))))) > 0)) %>% 
	mutate(PIB_final_in_pilot_mpra_ccre_variant = (PIB_final_in_pilot_adaptive_variant & PIB_final_mpra_ccre_variant))

mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	mutate(PIB_final_introgressed_variant = VariantID %in% introgressed_variants_tb$VariantID) %>% 
	mutate(PIB_final_adaptive_variant = VariantID %in% adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_final_corehaps_variant = VariantID %in% corehaps_variants_tb$VariantID) %>% 
	mutate(PIB_final_mpra_ccre_variant = VariantID %in% mpra_ccre_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_adaptive_variant = VariantID %in% pilot_adaptive_variants_tb$VariantID) %>% 
	mutate(PIB_pilot_mpra_ccre_variant = (VariantID %in% pilot_mpra_ccre_variants_tb$VariantID)) %>% 
	mutate(PIB_final_in_pilot_adaptive_variant = 
		(rowSums(!is.na(across(all_of(paste0("TractID_", map_pop_to_col_name[names(map_pop_to_pilot_pop)]))))) > 0)) %>% 
	mutate(PIB_final_in_pilot_mpra_ccre_variant = (PIB_final_in_pilot_adaptive_variant & PIB_final_mpra_ccre_variant))

# sort entries
adaptive_tracts_tb <- adaptive_tracts_tb %>% 
	GRanges() %>% sort() %>% as_tibble()
introgressed_tracts_tb <- introgressed_tracts_tb %>% 
	GRanges() %>% sort() %>% as_tibble()
introgressed_tracts_noOCN_tb <- introgressed_tracts_noOCN_tb %>% 
	GRanges() %>% sort() %>% as_tibble()
adaptive_corehaps_tb <- adaptive_corehaps_tb %>% 
	GRanges() %>% sort() %>% as_tibble()
introgressed_corehaps_tb <- introgressed_corehaps_tb %>% 
	GRanges() %>% sort() %>% as_tibble()
introgressed_variants_tb <- introgressed_variants_tb %>% 
	arrange(VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	GRanges() %>% sort() %>% as_tibble()
adaptive_variants_tb <- adaptive_variants_tb %>% 
	arrange(VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	GRanges() %>% sort() %>% as_tibble()
corehaps_variants_tb <- corehaps_variants_tb %>% 
	arrange(VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	GRanges() %>% sort() %>% as_tibble()
mpra_ccre_variants_tb <- mpra_ccre_variants_tb %>% 
	arrange(VariantCHROM, VariantPOS, VariantREF, VariantALT) %>% 
	GRanges() %>% sort() %>% as_tibble()

introgressed_nean_variants_tb <- introgressed_variants_tb %>% 
	filter(Alleles_NeandertalClass)
introgressed_deni_variants_tb <- introgressed_variants_tb %>% 
	filter(Alleles_DenisovanClass)
introgressed_ambig_variants_tb <- introgressed_variants_tb %>% 
	filter(Alleles_AmbiguousClass)

adaptive_nean_variants_tb <- adaptive_variants_tb %>% 
	filter(Alleles_NeandertalClass)
adaptive_deni_variants_tb <- adaptive_variants_tb %>% 
	filter(Alleles_DenisovanClass)
adaptive_ambig_variants_tb <- adaptive_variants_tb %>% 
	filter(Alleles_AmbiguousClass)

# create GRanges
adaptive_tracts_gr <- adaptive_tracts_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
introgressed_tracts_gr <- introgressed_tracts_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
introgressed_tracts_noOCN_gr <- introgressed_tracts_noOCN_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
adaptive_corehaps_gr <- adaptive_corehaps_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
introgressed_corehaps_gr <- introgressed_corehaps_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()

introgressed_variants_gr <- introgressed_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
adaptive_variants_gr <- adaptive_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
corehaps_variants_gr <- corehaps_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
mpra_ccre_variants_gr <- mpra_ccre_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()

introgressed_nean_variants_gr <- introgressed_nean_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
introgressed_deni_variants_gr <- introgressed_deni_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
introgressed_ambig_variants_gr <- introgressed_ambig_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()

adaptive_nean_variants_gr <- adaptive_nean_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
adaptive_deni_variants_gr <- adaptive_deni_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()
adaptive_ambig_variants_gr <- adaptive_ambig_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% sort()

introgressed_variants_lift38_gr <- introgressed_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
adaptive_variants_lift38_gr <- adaptive_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
corehaps_variants_lift38_gr <- corehaps_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()
mpra_ccre_variants_lift38_gr <- mpra_ccre_variants_tb %>% 
	GRanges() %>% seqinfo_fix("NCBI", "GRCh37") %>% 
	seqstyle_change("UCSC") %>% liftOver(hg19ToHg38) %>% unlist() %>% 
	seqinfo_fix("NCBI", "GRCh38") %>% sort()

# save as tables
write_tsv(adaptive_tracts_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_tracts.txt.gz"))
write_tsv(introgressed_tracts_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts.txt.gz"))
write_tsv(adaptive_corehaps_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_corehaps.txt.gz"))
write_tsv(introgressed_corehaps_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_corehaps.txt.gz"))
write_tsv(introgressed_variants_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.txt.gz"))
write_tsv(adaptive_variants_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants.txt.gz"))
write_tsv(corehaps_variants_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/corehaps_variants.txt.gz"))
write_tsv(mpra_ccre_variants_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants.txt.gz"))

# save as GRanges
saveRDS(adaptive_tracts_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_tracts.rds")
saveRDS(introgressed_tracts_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts.rds")
saveRDS(introgressed_tracts_noOCN_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts_noOCN.rds")
saveRDS(adaptive_corehaps_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_corehaps.rds")
saveRDS(introgressed_corehaps_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_corehaps.rds")
saveRDS(introgressed_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.rds")
saveRDS(adaptive_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants.rds")
saveRDS(corehaps_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/corehaps_variants.rds")
saveRDS(mpra_ccre_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants.rds")

# save as BED files
rtracklayer::export(adaptive_tracts_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_tracts.bed.gz", format="BED")
rtracklayer::export(introgressed_tracts_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts.bed.gz", format="BED")
rtracklayer::export(introgressed_tracts_noOCN_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_tracts_noOCN.bed.gz", format="BED")
rtracklayer::export(adaptive_corehaps_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_corehaps.bed.gz", format="BED")
rtracklayer::export(introgressed_corehaps_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_corehaps.bed.gz", format="BED")
rtracklayer::export(introgressed_variants_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants.bed.gz", format="BED")
rtracklayer::export(corehaps_variants_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/corehaps_variants.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants.bed.gz", format="BED")
rtracklayer::export(introgressed_variants_lift38_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/introgressed_variants_lift38.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_lift38_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/adaptive_variants_lift38.bed.gz", format="BED")
rtracklayer::export(corehaps_variants_lift38_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/corehaps_variants_lift38.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_lift38_gr %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/mpra_ccre_variants_lift38.bed.gz", format="BED")

# save as BED with chr
rtracklayer::export(adaptive_tracts_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/adaptive_tracts_chr.bed.gz", format="BED")
rtracklayer::export(introgressed_tracts_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/introgressed_tracts_chr.bed.gz", format="BED")
rtracklayer::export(introgressed_tracts_noOCN_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/introgressed_tracts_noOCN_chr.bed.gz", format="BED")
rtracklayer::export(adaptive_corehaps_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/adaptive_corehaps_chr.bed.gz", format="BED")
rtracklayer::export(introgressed_corehaps_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/introgressed_corehaps_chr.bed.gz", format="BED")
rtracklayer::export(introgressed_variants_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/introgressed_variants_chr.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/adaptive_variants_chr.bed.gz", format="BED")
rtracklayer::export(corehaps_variants_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/corehaps_variants_chr.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_ranges_chr/mpra_ccre_variants_chr.bed.gz", format="BED")
rtracklayer::export(introgressed_variants_lift38_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/introgressed_variants_lift38_chr.bed.gz", format="BED")
rtracklayer::export(adaptive_variants_lift38_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.bed.gz", format="BED")
rtracklayer::export(corehaps_variants_lift38_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/corehaps_variants_lift38_chr.bed.gz", format="BED")
rtracklayer::export(mpra_ccre_variants_lift38_gr %>% seqstyle_change("UCSC") %>% reduce() %>% sort(), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/mpra_ccre_variants_lift38_chr.bed.gz", format="BED")

# save as VCF files
convert_gr_to_vcf(introgressed_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants.vcf")
convert_gr_to_vcf(adaptive_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants.vcf")
convert_gr_to_vcf(corehaps_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/corehaps_variants.vcf")
convert_gr_to_vcf(mpra_ccre_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants.vcf")
convert_gr_to_vcf(introgressed_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(corehaps_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/corehaps_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/mpra_ccre_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(introgressed_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/introgressed_variants_lift38.vcf")
convert_gr_to_vcf(adaptive_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/adaptive_variants_lift38.vcf")
convert_gr_to_vcf(corehaps_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/corehaps_variants_lift38.vcf")
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/mpra_ccre_variants_lift38.vcf")
convert_gr_to_vcf(introgressed_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/introgressed_variants_lift38_clean.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/adaptive_variants_lift38_clean.vcf", clean=TRUE)
convert_gr_to_vcf(corehaps_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/corehaps_variants_lift38_clean.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38/mpra_ccre_variants_lift38_clean.vcf", clean=TRUE)

# save as VCF with chr
convert_gr_to_vcf(introgressed_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/introgressed_variants_chr.vcf")
convert_gr_to_vcf(adaptive_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/adaptive_variants_chr.vcf")
convert_gr_to_vcf(corehaps_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/corehaps_variants_chr.vcf")
convert_gr_to_vcf(mpra_ccre_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/mpra_ccre_variants_chr.vcf")
convert_gr_to_vcf(introgressed_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/introgressed_variants_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/adaptive_variants_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(corehaps_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/corehaps_variants_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_chr/mpra_ccre_variants_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(introgressed_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/introgressed_variants_lift38_chr.vcf")
convert_gr_to_vcf(adaptive_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/adaptive_variants_lift38_chr.vcf")
convert_gr_to_vcf(corehaps_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/corehaps_variants_lift38_chr.vcf")
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/mpra_ccre_variants_lift38_chr.vcf")
convert_gr_to_vcf(introgressed_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/introgressed_variants_lift38_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/adaptive_variants_lift38_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(corehaps_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/corehaps_variants_lift38_clean_chr.vcf", clean=TRUE)
convert_gr_to_vcf(mpra_ccre_variants_lift38_gr %>% seqstyle_change("UCSC"), "../../results/1a-preprocess_PIBv1_MPRA_final/alternate_variants_lift38_chr/mpra_ccre_variants_lift38_clean_chr.vcf", clean=TRUE)

# save flagged
introgressed_variants_flagged_liftOver_tb <- filter(introgressed_variants_tb, Flag_variant_failed_liftOver_hg38)
adaptive_variants_flagged_liftOver_tb <- filter(adaptive_variants_tb, Flag_variant_failed_liftOver_hg38)
corehaps_variants_flagged_liftOver_tb <- filter(corehaps_variants_tb, Flag_variant_failed_liftOver_hg38)
mpra_ccre_variants_flagged_liftOver_tb <- filter(mpra_ccre_variants_tb, Flag_variant_failed_liftOver_hg38)
write_tsv(introgressed_variants_flagged_liftOver_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/flagged_variants/introgressed_variants_flagged_liftOver.txt.gz"))
write_tsv(adaptive_variants_flagged_liftOver_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/flagged_variants/adaptive_variants_flagged_liftOver.txt.gz"))
write_tsv(corehaps_variants_flagged_liftOver_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/flagged_variants/corehaps_variants_flagged_liftOver.txt.gz"))
write_tsv(mpra_ccre_variants_flagged_liftOver_tb, gzfile("../../results/1a-preprocess_PIBv1_MPRA_final/flagged_variants/mpra_ccre_variants_flagged_liftOver.txt.gz"))

# save classes as VCFs
convert_gr_to_vcf(introgressed_nean_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_nean_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(introgressed_deni_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_deni_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(introgressed_ambig_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/introgressed_ambig_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_nean_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_nean_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_deni_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_deni_variants_clean.vcf", clean=TRUE)
convert_gr_to_vcf(adaptive_ambig_variants_gr, "../../results/1a-preprocess_PIBv1_MPRA_final/adaptive_ambig_variants_clean.vcf", clean=TRUE)
