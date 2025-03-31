#!/bin/R

# for data analysis
library(tidyverse)
library(data.table)
library(wrapr)

# for bio data analysis
library(plyranges)

find_overlap_map <- function(gr_query, gr_subject, names_query, names_subject, values_query="none", mode="none", mode_sep=";") {
	overlap <- findOverlaps(
		GRanges(gr_query), 
		GRanges(gr_subject)
	)
	let(c(VAR1=names_query, VAR2=names_subject), 
		mapping <- overlap %>% 
		as_tibble() %>% 
		mutate(
			VAR1 = gr_query$VAR1[queryHits(overlap)], 
			VAR2 = gr_subject$VAR2[subjectHits(overlap)]
		) %>% 
		dplyr::select(VAR1, VAR2)
	)
	if (mode=="left") {
		let(c(VAR1=names_query, VAR2=names_subject, VAR3=ifelse(values_query=="none", names_query, values_query)), 
			mapping <- mapping %>% 
				group_by(VAR2) %>% 
				dplyr::summarise(VAR3 = paste(VAR1, collapse=";"))
		)
	}
	if (mode=="right") {
		let(c(VAR1=names_query, VAR2=names_subject, VAR3=ifelse(values_query=="none", names_query, values_query)), 
			mapping <- mapping %>% 
				group_by(VAR2) %>% 
				dplyr::summarise(VAR3 = paste(VAR1, collapse=";"))
		)
	}
	return(mapping)
}

count_overlap <- function(gr_query, gr_subject, names_query, names_subject, values_query="none") {
	overlap <- findOverlaps(
		GRanges(gr_query), 
		GRanges(gr_subject)
	)
	let(c(VAR1=names_query, VAR2=names_subject), 
		mapping <- overlap %>% 
		as_tibble() %>% 
		mutate(
			VAR1 = gr_query$VAR1[queryHits(overlap)], 
			VAR2 = gr_subject$VAR2[subjectHits(overlap)]
		) %>% 
		dplyr::select(VAR1, VAR2)
	)
	let(c(VAR1=names_query, VAR2=names_subject, VAR3=ifelse(values_query=="none", names_query, values_query)), 
		mapping <- mapping %>% 
			group_by(VAR2) %>% 
			dplyr::summarise(VAR3 = dplyr::n())
	)
	return(mapping)
}
