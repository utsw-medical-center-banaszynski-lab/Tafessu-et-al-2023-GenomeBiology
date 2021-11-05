library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
library("rtracklayer")

WT <- readGAlignments("WT.bam")
H33KO <- readGAlignments("H33KO.bam")

WT.gr <- granges(WT)
H33KO.gr <- granges(H33KO)

##Assign to variable
WT <- WT.gr 
H33KO <- H33KO.gr

library_WT <- NROW(WT)
library_H33KO <- NROW(H33KO)

#Import BED files
bed_1 <- import("bed1.bed", format = "BED")

bed_1_rpkm_WT <- countOverlaps(bed_1, WT) / (width(bed_1)/1000 * library_WT/1000000)
bed_1_rpkm_H33KO <- countOverlaps(bed_1, H33KO) / (width(bed_1)/1000 * library_H33KO/1000000)

df_bed_1 <- data.frame(bed_1_rpkm_WT, bed_1_rpkm_H33KO)

pdf('box_plot.pdf')
boxplot(df_bed_1, col=(c("gray","red")), main="Title", ylab="RPKM",outline=FALSE, notch=TRUE, names=c("WT","H3.3KO"), yaxt="n", cex.axis=1,las=2,lwd=4,lty=1, cex.axis=1,las=2,lwd=4,lty=1, ylim=c(0,4))
axis(2, at=c(0,2,4))
dev.off()

#Wilcoxon rank sum test

wilcox.test(df_bed_1$bed_1_rpkm_WT, df_bed_1$bed_1_rpkm_H33KO,conf.int=TRUE)
