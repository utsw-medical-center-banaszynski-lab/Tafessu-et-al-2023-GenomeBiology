library("GenomicRanges")
library("GenomicAlignments")
library("rtracklayer")
library(readxl)
library(ggplot2)
library(ggrepel)

a6 <- readGAlignments("WT_GROseq.sorted.bam")
b6 <- readGAlignments("H33KO_GROseq.sorted.bam")


a6.gr <- granges(a6)
b6.gr <- granges(b6)


a6 <- a6.gr 
b6 <- b6.gr

library_a6 <- NROW(a6)
library_b6 <- NROW(b6)


#create list of bed file names, renamed to "motif.bed"

files <- list.files(path=".", pattern=".bed", all.files=F, full.names=F)
for (file in files) {
  motif <- import(file, format="BED")
}

for (file in files)  
{
#import motif BED, renamed to "motif.bed"
motif <- import(file, format="BED")
# Calculate RPKM
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_a6 <- countOverlaps(motif, a6) / (width(motif)/1000 * library_a6/1000000)
rpkm_b6 <- countOverlaps(motif, b6) / (width(motif)/1000 * library_b6/1000000)


rpkm_motif_1 <- data.frame(rpkm_a6, rpkm_b6)

#ratio
rpkm_motif_1[,"GRO_ratioKO.WT"] <- rpkm_motif_1[,2]/rpkm_motif_1[,1]

#log2ratio
rpkm_motif_1[,"GRO_log2ratioKO.WT"] <- log2(rpkm_motif_1[,3])


#filter out rows with "NA" or "Inf"
rpkm_motif_1 <- rpkm_motif_1 %>%
  filter_if(is.numeric, all_vars(!is.na(.))) %>%
  filter_if(is.numeric, all_vars(!is.infinite(.)))

#add median to "FC_GRO.txt"
med_FC <- median (rpkm_motif_1[,3])
med_log2FC <- median (rpkm_motif_1[,4])
write(med_log2FC,file="FC_GRO.txt",append=TRUE)
write.table(rpkm_motif_1, file=paste(myfunc(file), "FC_GRO.xls", sep = "_"), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

}


###########################################################################################################

#Make scatterplot from "TF_GRO_vs_DBS.xls" containing columns for TF name, GRO-seq log2FC, Differential Binding Score (DBS) and label ("yes" or "no")

Joint_TF_label <- read.table(file = "TF_GRO_vs_DBS.xls", header= T, sep="\t",stringsAsFactors=FALSE)
head(Joint_TF_label)


corr_plot <- ggplot(Joint_TF_label, aes(x=DBS, y=log2FC_GRO, color=ToLabel)) + geom_point()+ scale_color_manual(values = c("no" = "gray", "yes" = "forestgreen")) + xlim(-0.5,0) + ylim(-0.6, -0.2) + theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


corr_plot <- corr_plot + geom_label_repel(data = subset(Joint_TF_label, ToLabel == "yes"), 
                                          aes(label = TF_ID, size = NULL, color = NULL),
                                          nudge_y = -0.1,
                                          segment.size  = 0.2,
                                          segment.color = "grey50",
                                          direction     = "x"
)

corr_plot_fit <- corr_plot + geom_smooth(method=lm,linetype=2,colour="red",se=F)

pdf("scatterplot.pdf")
corr_plot
dev.off()
