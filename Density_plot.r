library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")

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


