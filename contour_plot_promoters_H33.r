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
h33k <- readGAlignments("H33.bam")

a1.gr <- granges(a1)
b1.gr <- granges(b1)
a2.gr <- granges(a2)
b2.gr <- granges(b2)
a3.gr <- granges(a3)
b3.gr <- granges(b3)
a4.gr <- granges(a4)
b4.gr <- granges(b4)
h33k.gr <- granges(h33k)


a1 <- a1.gr 
b1 <- b1.gr
a2 <- a2.gr 
b2 <- b2.gr
a3 <- a3.gr 
b3 <- b3.gr
a4 <- a4.gr 
b4 <- b4.gr
h33k <- h33k.gr

library_a1 <- NROW(a1)
library_b1 <- NROW(b1)
library_a2 <- NROW(a2)
library_b2 <- NROW(b2)
library_a3 <- NROW(a3)
library_b3 <- NROW(b3)
library_a4 <- NROW(a4)
library_b4 <- NROW(b4)
library_h33k <- NROW(h33k)


promoters <- import("mESC_expressed_gene_promoters.bed", format = "BED")

###########################################################################################################

#p300

# Calculate RPKM
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a1 <- countOverlaps(promoters, a1) / (width(promoters)/1000 * library_a1/1000000)
rpkm_b1 <- countOverlaps(promoters, b1) / (width(promoters)/1000 * library_b1/1000000)

rpkm_h33k <- countOverlaps(promoters, h33k) / (width(promoters)/1000 * library_h33k/1000000)

rpkm_promoters_1 <- data.frame(rpkm_a1, rpkm_b1, rpkm_h33k)

#ratio
rpkm_promoters_1[,"p300_ratioKO.WT"] <- rpkm_promoters_1[,2]/rpkm_promoters_1[,1]

#log2ratio
rpkm_promoters_1[,"p300_log2ratioKO.WT"] <- log2(rpkm_promoters_1[,4])

pdf("contourplot_h33k_log2p300.pdf")
ggplot(data=rpkm_promoters_1,aes(rpkm_h33k,p300_log2ratioKO.WT)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(0, 2) + ylim(-1, 1) + theme(aspect.ratio=1)
guides(alpha="none")
dev.off()


###########################################################################################################

#H3K4me3

# Calculate RPKM
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a2 <- countOverlaps(promoters, a2) / (width(promoters)/1000 * library_a2/1000000)
rpkm_b2 <- countOverlaps(promoters, b2) / (width(promoters)/1000 * library_b2/1000000)

rpkm_h33k <- countOverlaps(promoters, h33k) / (width(promoters)/1000 * library_h33k/1000000)

rpkm_promoters_2 <- data.frame(rpkm_a2, rpkm_b2, rpkm_h33k)

#ratio
rpkm_promoters_2[,"H3K4me3_ratioKO.WT"] <- rpkm_promoters_2[,2]/rpkm_promoters_2[,1]

#log2ratio
rpkm_promoters_2[,"H3K4me3_log2ratioKO.WT"] <- log2(rpkm_promoters_2[,4])

pdf("contourplot_h33k_log2H3K4me3.pdf")
ggplot(data=rpkm_promoters_2,aes(rpkm_h33k,H3K4me3_log2ratioKO.WT)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(0, 2) + ylim(-1, 1) + theme(aspect.ratio=1)
guides(alpha="none")
dev.off()

###########################################################################################################

#H3K27ac

# Calculate RPKM
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a3 <- countOverlaps(promoters, a3) / (width(promoters)/1000 * library_a3/1000000)
rpkm_b3 <- countOverlaps(promoters, b3) / (width(promoters)/1000 * library_b3/1000000)

rpkm_h33k <- countOverlaps(promoters, h33k) / (width(promoters)/1000 * library_h33k/1000000)

rpkm_promoters_3 <- data.frame(rpkm_a3, rpkm_b3, rpkm_h33k)

#ratio
rpkm_promoters_3[,"H3K27ac_ratioKO.WT"] <- rpkm_promoters_3[,2]/rpkm_promoters_3[,1]

#log2ratio
rpkm_promoters_3[,"H3K27ac_log2ratioKO.WT"] <- log2(rpkm_promoters_3[,4])

pdf("contourplot_h33k_log2H3K27ac.pdf")
ggplot(data=rpkm_promoters_3,aes(rpkm_h33k,H3K27ac_log2ratioKO.WT)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(0, 2) + ylim(-2, 2) + theme(aspect.ratio=1)
guides(alpha="none")
dev.off()

###########################################################################################################

#BRD4

# Calculate RPKM
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a4 <- countOverlaps(promoters, a4) / (width(promoters)/1000 * library_a4/1000000)
rpkm_b4 <- countOverlaps(promoters, b4) / (width(promoters)/1000 * library_b4/1000000)

rpkm_h33k <- countOverlaps(promoters, h33k) / (width(promoters)/1000 * library_h33k/1000000)

rpkm_promoters_4 <- data.frame(rpkm_a4, rpkm_b4, rpkm_h33k)


#ratio
rpkm_promoters_4[,"BRD4_ratioKO.WT"] <- rpkm_promoters_4[,2]/rpkm_promoters_4[,1]

#log2ratio
rpkm_promoters_4[,"BRD4_log2ratioKO.WT"] <- log2(rpkm_promoters_4[,4])

pdf("contourplot_h33k_log2BRD4.pdf")
ggplot(data=rpkm_promoters_4,aes(rpkm_h33k,BRD4_log2ratioKO.WT)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(0, 2) + ylim(-2, 2) + theme(aspect.ratio=1)
guides(alpha="none")
dev.off()

###########################################################################################################

# H33 vs ATAC in WT

rpkm_promoters_6 <- data.frame(rpkm_a5, rpkm_h33k)

pdf("contourplot_WT_h33k_ATAC.pdf")
ggplot(data=rpkm_promoters_6,aes(rpkm_h33k,rpkm_a5)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..), geom='polygon',colour='black') + 
  scale_fill_continuous(low="green",high="red") +
  geom_smooth(method=lm,linetype=2,colour="red",se=F) + xlim(0, 2) + ylim(0, 3) + theme(aspect.ratio=1)
guides(alpha="none")
dev.off()

