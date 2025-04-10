#!/bin/R
library(tidyverse)
a <- read_tsv("../results/MPRAmodel_Tucci_MPRA_pilot_K562_K562_emVAR_glm_20230106.out")
b <- a %>% filter(Skew_logFDR_act > -log10(0.05))
pdf("MPRAmodel_Tucci_MPRA_pilot_K562_emVar.pdf")
plot(pmax((a$A_logPadj_BH), (a$B_logPadj_BH)), a$Log2Skew, xlab="K562 expression -log10 FDR (max allele)", ylab="K562 allelic log2 Skew (ALT/REF)")
points(pmax((b$A_logPadj_BH), (b$B_logPadj_BH)), b$Log2Skew, col="#367db6")
abline(h = 0, col="#eeac6f")
legend('topright',c('Not Significant', 'emVar'),   
	col=c("black", "#367db6"), horiz=F, cex=0.8,
	pch=c(20,20), bg='white')
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
dev.off()

pdf("MPRAmodel_Tucci_MPRA_pilot_K562_emVar_2.pdf")
plot(a$Log2Skew, a$Skew_logFDR_act, ylab="K562 allelic skew -log10 FDR", xlab="K562 allelic log2 Skew (ALT/REF)")
points(b$Log2Skew, b$Skew_logFDR_act, col="#367db6")
# abline(h = 0, col="#eeac6f")
abline(v = 0, col="#eeac6f")
legend('topright',c('Not Significant', 'emVar'),   
	col=c("black", "#367db6"), horiz=F, cex=0.8,
	pch=c(20,20), bg='white')
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
dev.off()

b_ns <- a %>% filter(Skew_logFDR_act < -log10(0.05))
pdf("MPRAmodel_Tucci_MPRA_pilot_K562_emVar_3.pdf")
plot(a$Log2Skew, a$Skew_logFDR_act, ylab="K562 allelic skew -log10 FDR", xlab="K562 allelic log2 Skew (ALT/REF)", col="#367db6")
points(b_ns$Log2Skew, b_ns$Skew_logFDR_act, col="black")
# abline(h = 0, col="#eeac6f")
abline(v = 0, col="#eeac6f")
legend('topright',c('Not Significant', 'emVar'),   
	col=c("black", "#367db6"), horiz=F, cex=0.8,
	pch=c(20,20), bg='white')
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
dev.off()
