library("GenomicRanges")
library("GenomicAlignments")
library("rtracklayer")
library("ggplot2")

#read BAM files for ChIP data to be compared between WT and H3.3 KO, and H3.3 ChIP in WT.
a1 <- readGAlignments("WT.bam")
b1 <- readGAlignments("H33KO.bam")
h33 <- readGAlignments("H33_ChIP.bam")

a1.gr <- granges(a1)
b1.gr <- granges(b1)
h33.gr <- granges(h33)

a1 <- a1.gr 
b1 <- b1.gr
h33 <- h33.gr

library_a1 <- NROW(a1)
library_b1 <- NROW(b1)
library_h33 <- NROW(h33)

promoters <- import("mESC_expressed_gene_promoters.bed", format = "BED")

###########################################################################################################

# Calculate RPKM
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a1 <- countOverlaps(promoters, a1) / (width(promoters)/1000 * library_a1/1000000)
rpkm_b1 <- countOverlaps(promoters, b1) / (width(promoters)/1000 * library_b1/1000000)

rpkm_h33 <- countOverlaps(promoters, h33) / (width(promoters)/1000 * library_h33/1000000)

rpkm_promoters_1 <- data.frame(rpkm_a1, rpkm_b1, rpkm_h33)

#ratio
rpkm_promoters_1[,"ratioKO.WT"] <- rpkm_promoters_1[,2]/rpkm_promoters_1[,1]

#log2ratio
rpkm_promoters_1[,"log2ratioKO.WT"] <- log2(rpkm_promoters_1[,4])

pdf("contourplot.pdf")
ggplot(data=rpkm_promoters_1,aes(rpkm_h33,p300_log2ratioKO.WT)) + stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + scale_fill_continuous(low="green",high="red") + geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(0, 2) + ylim(-1, 1) + theme(aspect.ratio=1)
guides(alpha="none")
dev.off()
