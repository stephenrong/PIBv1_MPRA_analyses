#Read in Count and Attributes tables
source("MPRAmodel_edited.R")

count_PROJ <- read.delim("../MPRAcount_output_Jurkat/Tucci_MPRA_pilot_Jurkat_Jurkat.counts", stringsAsFactors=F)
attr_PROJ <- read.delim("MPRAmodel_Tucci_MPRA_pilot_attributes.txt", stringsAsFactors=F)
cond_PROJ <- read.delim("MPRAmodel_Tucci_MPRA_pilot_condition_Jurkat.txt", stringsAsFactors=F, row.names=1)  ###

MPRAmodel(count_PROJ, attr_PROJ, cond_PROJ, filePrefix="MPRAmodel_Tucci_MPRA_pilot_Jurkat", projectName="MPRAmodel_Tucci_MPRA_pilot_Jurkat", negCtrlName="Tucci_MPRA_ctrl_neg", posCtrlName="Tucci_MPRA_ctrl_pos", emvarCtrlName="Tucci_MPRA_ctrl_emvar")  ###

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
