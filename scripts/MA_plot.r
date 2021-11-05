library(gdata)
library(gplots)
library(reshape)
library(R.utils)
library(dplyr)
library(readxl)
library(ggplot2)
library(plyr)
library(devtools)
source_gist("524eade46135f6348140")
library(ggpubr)
library(ggrepel)

TF_list <- read.table(file = "TF_cofactor_list.xls", header= T, sep="\t",stringsAsFactors=FALSE)

ma_plot <- ggplot(TF_list, aes(x=log2baseMean, y=log2FoldChange, color=Type)) + geom_point(size=0.5)+ scale_color_manual(values = c("no" = "gray", "TF" = "forestgreen", "cofactor" = "purple")) + xlim(0,16) + ylim(-10, 10) + theme(aspect.ratio=1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black") + scale_y_continuous(breaks = seq(-10, 10, by = 2)))

ma_plot <- ma_plot + geom_label_repel(data = subset(TF_list, Type == "cofactor"), 
                             aes(label = Symbol, size = NULL, color = NULL),
                             nudge_y = -2,
                             segment.size  = 0.2,
                             segment.color = "grey50",
                             direction     = "x"
)

ma_plot <- ma_plot + geom_label_repel(data = subset(TF_list, Type == "TF"), 
                                      aes(label = Symbol, size = NULL, color = NULL),
                                      nudge_y = 2,
                                      segment.size  = 0.2,
                                      segment.color = "grey50",
                                      direction     = "x"
)

pdf("MAplot.pdf")
ma_plot
dev.off()
