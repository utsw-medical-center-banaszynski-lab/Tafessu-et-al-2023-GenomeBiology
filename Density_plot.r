#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")
#options(mc.cores=getCores(4))
#install.packages("vioplot")
#library(vioplot)


a1 <- readGAlignments("WT_p300.bam")
b1 <- readGAlignments("H33KO_p300.bam")
a2 <- readGAlignments("WT_H3K4me3.bam")
b2 <- readGAlignments("H33KO_H3K4me3.bam")
a3 <- readGAlignments("WT_H3K27ac.bam")
b3 <- readGAlignments("H33KO_H3K27ac.bam")
a4 <- readGAlignments("WT_BRD4.bam")
b4 <- readGAlignments("H33KO_BRD4.bam")


a1.gr <- granges(a1)
b1.gr <- granges(b1)
a2.gr <- granges(a2)
b2.gr <- granges(b2)
a3.gr <- granges(a3)
b3.gr <- granges(b3)
a4.gr <- granges(a4)
b4.gr <- granges(b4)

##Assign to variable
a1 <- a1.gr 
b1 <- b1.gr
a2 <- a2.gr 
b2 <- b2.gr
a3 <- a3.gr 
b3 <- b3.gr
a4 <- a4.gr 
b4 <- b4.gr


library_a1 <- NROW(a1)
library_b1 <- NROW(b1)
library_a2 <- NROW(a2)
library_b2 <- NROW(b2)
library_a3 <- NROW(a3)
library_b3 <- NROW(b3)
library_a4 <- NROW(a4)
library_b4 <- NROW(b4)

KLF4 <- import("KLF4_bound_promoters.bed", format = "BED")
KLF5 <- import("KLF5_bound_promoters.bed", format = "BED")
POU5F1 <- import("POU5F1_bound_promoters.bed", format = "BED")
STAT3 <- import("STAT3_bound_promoters.bed", format = "BED")
ZFP42 <- import("ZFP42_bound_promoters.bed", format = "BED")


########################################################################################################
#KLF4
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
KLF4_a1 <- countOverlaps(KLF4, a1) / (width(KLF4)/1000 * library_a1/1000000)
KLF4_b1 <- countOverlaps(KLF4, b1) / (width(KLF4)/1000 * library_b1/1000000)
KLF4_a2 <- countOverlaps(KLF4, a2) / (width(KLF4)/1000 * library_a2/1000000)
KLF4_b2 <- countOverlaps(KLF4, b2) / (width(KLF4)/1000 * library_b2/1000000)
KLF4_a3 <- countOverlaps(KLF4, a3) / (width(KLF4)/1000 * library_a3/1000000)
KLF4_b3 <- countOverlaps(KLF4, b3) / (width(KLF4)/1000 * library_b3/1000000)
KLF4_a4 <- countOverlaps(KLF4, a4) / (width(KLF4)/1000 * library_a4/1000000)
KLF4_b4 <- countOverlaps(KLF4, b4) / (width(KLF4)/1000 * library_b4/1000000)


rpkm_KLF4_p300 <- data.frame(KLF4_a1, KLF4_b1)
rpkm_KLF4_p300[,"log2ratioKO.WT"] <- log2(rpkm_KLF4_p300[,2]/rpkm_KLF4_p300[,1])

rpkm_KLF4_H3K4me3 <- data.frame(KLF4_a2, KLF4_b2)
rpkm_KLF4_H3K4me3[,"log2ratioKO.WT"] <- log2(rpkm_KLF4_H3K4me3[,2]/rpkm_KLF4_H3K4me3[,1])

rpkm_KLF4_H3K27ac <- data.frame(KLF4_a3, KLF4_b3)
rpkm_KLF4_H3K27ac[,"log2ratioKO.WT"] <- log2(rpkm_KLF4_H3K27ac[,2]/rpkm_KLF4_H3K27ac[,1])

rpkm_KLF4_BRD4 <- data.frame(KLF4_a4, KLF4_b4)
rpkm_KLF4_BRD4[,"log2ratioKO.WT"] <- log2(rpkm_KLF4_BRD4[,2]/rpkm_KLF4_BRD4[,1])

pdf("densityplot_H33WT_H33KO_KLF4_bound_promoters.pdf")
plot(density(rpkm_KLF4_H3K4me3$log2ratioKO.WT ,na.rm=TRUE),col="darkorange", lty=1, lwd=4, main="KLF4-bound active promoters (H3.3 KO/WT)", xlab="log2 RPKM KO/WT", ylim=c(0,2.5), xlim=c(-3,3))
lines(density(rpkm_KLF4_p300$log2ratioKO.WT ,na.rm=TRUE),col="olivedrab", lty=1, lwd=4)
lines(density(rpkm_KLF4_H3K27ac$log2ratioKO.WT ,na.rm=TRUE),col="dodgerblue3", lty=1, lwd=4)
lines(density(rpkm_KLF4_BRD4$log2ratioKO.WT ,na.rm=TRUE),col="darkorchid4", lty=1, lwd=4)
legend('topright',c('H3K4me3','p300', 'H3K27ac', 'BRD4'), fill = c("darkorange","olivedrab","dodgerblue3","darkorchid4"), bty = 'n', border = NA)
abline(v = mean(rpkm_KLF4_H3K4me3$log2ratioKO.WT), col="darkorange", lwd=2, lty=2)
abline(v = mean(rpkm_KLF4_p300$log2ratioKO.WT), col="olivedrab", lwd=2, lty=2)
abline(v = mean(rpkm_KLF4_H3K27ac$log2ratioKO.WT), col="dodgerblue3", lwd=2, lty=2)
abline(v = mean(rpkm_KLF4_BRD4$log2ratioKO.WT), col="darkorchid4", lwd=2, lty=2)
dev.off()
########################################################################################################

#KLF5
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
KLF5_a1 <- countOverlaps(KLF5, a1) / (width(KLF5)/1000 * library_a1/1000000)
KLF5_b1 <- countOverlaps(KLF5, b1) / (width(KLF5)/1000 * library_b1/1000000)
KLF5_a2 <- countOverlaps(KLF5, a2) / (width(KLF5)/1000 * library_a2/1000000)
KLF5_b2 <- countOverlaps(KLF5, b2) / (width(KLF5)/1000 * library_b2/1000000)
KLF5_a3 <- countOverlaps(KLF5, a3) / (width(KLF5)/1000 * library_a3/1000000)
KLF5_b3 <- countOverlaps(KLF5, b3) / (width(KLF5)/1000 * library_b3/1000000)
KLF5_a4 <- countOverlaps(KLF5, a4) / (width(KLF5)/1000 * library_a4/1000000)
KLF5_b4 <- countOverlaps(KLF5, b4) / (width(KLF5)/1000 * library_b4/1000000)


rpkm_KLF5_p300 <- data.frame(KLF5_a1, KLF5_b1)
rpkm_KLF5_p300[,"log2ratioKO.WT"] <- log2(rpkm_KLF5_p300[,2]/rpkm_KLF5_p300[,1])

rpkm_KLF5_H3K4me3 <- data.frame(KLF5_a2, KLF5_b2)
rpkm_KLF5_H3K4me3[,"log2ratioKO.WT"] <- log2(rpkm_KLF5_H3K4me3[,2]/rpkm_KLF5_H3K4me3[,1])

rpkm_KLF5_H3K27ac <- data.frame(KLF5_a3, KLF5_b3)
rpkm_KLF5_H3K27ac[,"log2ratioKO.WT"] <- log2(rpkm_KLF5_H3K27ac[,2]/rpkm_KLF5_H3K27ac[,1])

rpkm_KLF5_BRD4 <- data.frame(KLF5_a4, KLF5_b4)
rpkm_KLF5_BRD4[,"log2ratioKO.WT"] <- log2(rpkm_KLF5_BRD4[,2]/rpkm_KLF5_BRD4[,1])

pdf("densityplot_H33WT_H33KO_KLF5_bound_promoters.pdf")
plot(density(rpkm_KLF5_H3K4me3$log2ratioKO.WT ,na.rm=TRUE),col="darkorange", lty=1, lwd=4, main="KLF5-bound active promoters (H3.3 KO/WT)", xlab="log2 RPKM KO/WT", ylim=c(0,2.5), xlim=c(-3,3))
lines(density(rpkm_KLF5_p300$log2ratioKO.WT ,na.rm=TRUE),col="olivedrab", lty=1, lwd=4)
lines(density(rpkm_KLF5_H3K27ac$log2ratioKO.WT ,na.rm=TRUE),col="dodgerblue3", lty=1, lwd=4)
lines(density(rpkm_KLF5_BRD4$log2ratioKO.WT ,na.rm=TRUE),col="darkorchid4", lty=1, lwd=4)
legend('topright',c('H3K4me3','p300', 'H3K27ac', 'BRD4'), fill = c("darkorange","olivedrab","dodgerblue3","darkorchid4"), bty = 'n', border = NA)
abline(v = mean(rpkm_KLF5_H3K4me3$log2ratioKO.WT), col="darkorange", lwd=2, lty=2)
abline(v = mean(rpkm_KLF5_p300$log2ratioKO.WT), col="olivedrab", lwd=2, lty=2)
abline(v = mean(rpkm_KLF5_H3K27ac$log2ratioKO.WT), col="dodgerblue3", lwd=2, lty=2)
abline(v = mean(rpkm_KLF5_BRD4$log2ratioKO.WT), col="darkorchid4", lwd=2, lty=2)
dev.off()
########################################################################################################

#POU5F1
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
POU5F1_a1 <- countOverlaps(POU5F1, a1) / (width(POU5F1)/1000 * library_a1/1000000)
POU5F1_b1 <- countOverlaps(POU5F1, b1) / (width(POU5F1)/1000 * library_b1/1000000)
POU5F1_a2 <- countOverlaps(POU5F1, a2) / (width(POU5F1)/1000 * library_a2/1000000)
POU5F1_b2 <- countOverlaps(POU5F1, b2) / (width(POU5F1)/1000 * library_b2/1000000)
POU5F1_a3 <- countOverlaps(POU5F1, a3) / (width(POU5F1)/1000 * library_a3/1000000)
POU5F1_b3 <- countOverlaps(POU5F1, b3) / (width(POU5F1)/1000 * library_b3/1000000)
POU5F1_a4 <- countOverlaps(POU5F1, a4) / (width(POU5F1)/1000 * library_a4/1000000)
POU5F1_b4 <- countOverlaps(POU5F1, b4) / (width(POU5F1)/1000 * library_b4/1000000)


rpkm_POU5F1_p300 <- data.frame(POU5F1_a1, POU5F1_b1)
rpkm_POU5F1_p300[,"log2ratioKO.WT"] <- log2(rpkm_POU5F1_p300[,2]/rpkm_POU5F1_p300[,1])

rpkm_POU5F1_H3K4me3 <- data.frame(POU5F1_a2, POU5F1_b2)
rpkm_POU5F1_H3K4me3[,"log2ratioKO.WT"] <- log2(rpkm_POU5F1_H3K4me3[,2]/rpkm_POU5F1_H3K4me3[,1])

rpkm_POU5F1_H3K27ac <- data.frame(POU5F1_a3, POU5F1_b3)
rpkm_POU5F1_H3K27ac[,"log2ratioKO.WT"] <- log2(rpkm_POU5F1_H3K27ac[,2]/rpkm_POU5F1_H3K27ac[,1])

rpkm_POU5F1_BRD4 <- data.frame(POU5F1_a4, POU5F1_b4)
rpkm_POU5F1_BRD4[,"log2ratioKO.WT"] <- log2(rpkm_POU5F1_BRD4[,2]/rpkm_POU5F1_BRD4[,1])

pdf("densityplot_H33WT_H33KO_POU5F1_bound_promoters.pdf")
plot(density(rpkm_POU5F1_H3K4me3$log2ratioKO.WT ,na.rm=TRUE),col="darkorange", lty=1, lwd=4, main="POU5F1-bound active promoters (H3.3 KO/WT)", xlab="log2 RPKM KO/WT", ylim=c(0,2.5), xlim=c(-3,3))
lines(density(rpkm_POU5F1_p300$log2ratioKO.WT ,na.rm=TRUE),col="olivedrab", lty=1, lwd=4)
lines(density(rpkm_POU5F1_H3K27ac$log2ratioKO.WT ,na.rm=TRUE),col="dodgerblue3", lty=1, lwd=4)
lines(density(rpkm_POU5F1_BRD4$log2ratioKO.WT ,na.rm=TRUE),col="darkorchid4", lty=1, lwd=4)
legend('topright',c('H3K4me3','p300', 'H3K27ac', 'BRD4'), fill = c("darkorange","olivedrab","dodgerblue3","darkorchid4"), bty = 'n', border = NA)
abline(v = mean(rpkm_POU5F1_H3K4me3$log2ratioKO.WT), col="darkorange", lwd=2, lty=2)
abline(v = mean(rpkm_POU5F1_p300$log2ratioKO.WT), col="olivedrab", lwd=2, lty=2)
abline(v = mean(rpkm_POU5F1_H3K27ac$log2ratioKO.WT), col="dodgerblue3", lwd=2, lty=2)
abline(v = mean(rpkm_POU5F1_BRD4$log2ratioKO.WT), col="darkorchid4", lwd=2, lty=2)
dev.off()
########################################################################################################

#STAT3
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
STAT3_a1 <- countOverlaps(STAT3, a1) / (width(STAT3)/1000 * library_a1/1000000)
STAT3_b1 <- countOverlaps(STAT3, b1) / (width(STAT3)/1000 * library_b1/1000000)
STAT3_a2 <- countOverlaps(STAT3, a2) / (width(STAT3)/1000 * library_a2/1000000)
STAT3_b2 <- countOverlaps(STAT3, b2) / (width(STAT3)/1000 * library_b2/1000000)
STAT3_a3 <- countOverlaps(STAT3, a3) / (width(STAT3)/1000 * library_a3/1000000)
STAT3_b3 <- countOverlaps(STAT3, b3) / (width(STAT3)/1000 * library_b3/1000000)
STAT3_a4 <- countOverlaps(STAT3, a4) / (width(STAT3)/1000 * library_a4/1000000)
STAT3_b4 <- countOverlaps(STAT3, b4) / (width(STAT3)/1000 * library_b4/1000000)


rpkm_STAT3_p300 <- data.frame(STAT3_a1, STAT3_b1)
rpkm_STAT3_p300[,"log2ratioKO.WT"] <- log2(rpkm_STAT3_p300[,2]/rpkm_STAT3_p300[,1])

rpkm_STAT3_H3K4me3 <- data.frame(STAT3_a2, STAT3_b2)
rpkm_STAT3_H3K4me3[,"log2ratioKO.WT"] <- log2(rpkm_STAT3_H3K4me3[,2]/rpkm_STAT3_H3K4me3[,1])

rpkm_STAT3_H3K27ac <- data.frame(STAT3_a3, STAT3_b3)
rpkm_STAT3_H3K27ac[,"log2ratioKO.WT"] <- log2(rpkm_STAT3_H3K27ac[,2]/rpkm_STAT3_H3K27ac[,1])

rpkm_STAT3_BRD4 <- data.frame(STAT3_a4, STAT3_b4)
rpkm_STAT3_BRD4[,"log2ratioKO.WT"] <- log2(rpkm_STAT3_BRD4[,2]/rpkm_STAT3_BRD4[,1])

pdf("densityplot_H33WT_H33KO_STAT3_bound_promoters.pdf")
plot(density(rpkm_STAT3_H3K4me3$log2ratioKO.WT ,na.rm=TRUE),col="darkorange", lty=1, lwd=4, main="STAT3-bound active promoters (H3.3 KO/WT)", xlab="log2 RPKM KO/WT", ylim=c(0,2.5), xlim=c(-3,3))
lines(density(rpkm_STAT3_p300$log2ratioKO.WT ,na.rm=TRUE),col="olivedrab", lty=1, lwd=4)
lines(density(rpkm_STAT3_H3K27ac$log2ratioKO.WT ,na.rm=TRUE),col="dodgerblue3", lty=1, lwd=4)
lines(density(rpkm_STAT3_BRD4$log2ratioKO.WT ,na.rm=TRUE),col="darkorchid4", lty=1, lwd=4)
legend('topright',c('H3K4me3','p300', 'H3K27ac', 'BRD4'), fill = c("darkorange","olivedrab","dodgerblue3","darkorchid4"), bty = 'n', border = NA)
abline(v = mean(rpkm_STAT3_H3K4me3$log2ratioKO.WT), col="darkorange", lwd=2, lty=2)
abline(v = mean(rpkm_STAT3_p300$log2ratioKO.WT), col="olivedrab", lwd=2, lty=2)
abline(v = mean(rpkm_STAT3_H3K27ac$log2ratioKO.WT), col="dodgerblue3", lwd=2, lty=2)
abline(v = mean(rpkm_STAT3_BRD4$log2ratioKO.WT), col="darkorchid4", lwd=2, lty=2)
dev.off()
########################################################################################################

#ZFP42
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
ZFP42_a1 <- countOverlaps(ZFP42, a1) / (width(ZFP42)/1000 * library_a1/1000000)
ZFP42_b1 <- countOverlaps(ZFP42, b1) / (width(ZFP42)/1000 * library_b1/1000000)
ZFP42_a2 <- countOverlaps(ZFP42, a2) / (width(ZFP42)/1000 * library_a2/1000000)
ZFP42_b2 <- countOverlaps(ZFP42, b2) / (width(ZFP42)/1000 * library_b2/1000000)
ZFP42_a3 <- countOverlaps(ZFP42, a3) / (width(ZFP42)/1000 * library_a3/1000000)
ZFP42_b3 <- countOverlaps(ZFP42, b3) / (width(ZFP42)/1000 * library_b3/1000000)
ZFP42_a4 <- countOverlaps(ZFP42, a4) / (width(ZFP42)/1000 * library_a4/1000000)
ZFP42_b4 <- countOverlaps(ZFP42, b4) / (width(ZFP42)/1000 * library_b4/1000000)


rpkm_ZFP42_p300 <- data.frame(ZFP42_a1, ZFP42_b1)
rpkm_ZFP42_p300[,"log2ratioKO.WT"] <- log2(rpkm_ZFP42_p300[,2]/rpkm_ZFP42_p300[,1])

rpkm_ZFP42_H3K4me3 <- data.frame(ZFP42_a2, ZFP42_b2)
rpkm_ZFP42_H3K4me3[,"log2ratioKO.WT"] <- log2(rpkm_ZFP42_H3K4me3[,2]/rpkm_ZFP42_H3K4me3[,1])

rpkm_ZFP42_H3K27ac <- data.frame(ZFP42_a3, ZFP42_b3)
rpkm_ZFP42_H3K27ac[,"log2ratioKO.WT"] <- log2(rpkm_ZFP42_H3K27ac[,2]/rpkm_ZFP42_H3K27ac[,1])

rpkm_ZFP42_BRD4 <- data.frame(ZFP42_a4, ZFP42_b4)
rpkm_ZFP42_BRD4[,"log2ratioKO.WT"] <- log2(rpkm_ZFP42_BRD4[,2]/rpkm_ZFP42_BRD4[,1])

pdf("densityplot_H33WT_H33KO_ZFP42_bound_promoters.pdf")
plot(density(rpkm_ZFP42_H3K4me3$log2ratioKO.WT ,na.rm=TRUE),col="darkorange", lty=1, lwd=4, main="ZFP42-bound active promoters (H3.3 KO/WT)", xlab="log2 RPKM KO/WT", ylim=c(0,2.5), xlim=c(-3,3))
lines(density(rpkm_ZFP42_p300$log2ratioKO.WT ,na.rm=TRUE),col="olivedrab", lty=1, lwd=4)
lines(density(rpkm_ZFP42_H3K27ac$log2ratioKO.WT ,na.rm=TRUE),col="dodgerblue3", lty=1, lwd=4)
lines(density(rpkm_ZFP42_BRD4$log2ratioKO.WT ,na.rm=TRUE),col="darkorchid4", lty=1, lwd=4)
legend('topright',c('H3K4me3','p300', 'H3K27ac', 'BRD4'), fill = c("darkorange","olivedrab","dodgerblue3","darkorchid4"), bty = 'n', border = NA)
abline(v = mean(rpkm_ZFP42_H3K4me3$log2ratioKO.WT), col="darkorange", lwd=2, lty=2)
abline(v = mean(rpkm_ZFP42_p300$log2ratioKO.WT), col="olivedrab", lwd=2, lty=2)
abline(v = mean(rpkm_ZFP42_H3K27ac$log2ratioKO.WT), col="dodgerblue3", lwd=2, lty=2)
abline(v = mean(rpkm_ZFP42_BRD4$log2ratioKO.WT), col="darkorchid4", lwd=2, lty=2)
dev.off()
########################################################################################################
