require(rtracklayer)
require(ChIPpeakAnno)
require(Vennerable)



#import bed files

WT <- import.bed("WT.bed")
KO <- import.bed("KO.bed")

# compute overlaps
overlap <- findOverlaps(WT, KO)

# get subsets of binding sites
subset <- subsetByOverlaps(WT, KO)

# number of overlaps
length(subset)

#Make a Venn-diagram

WT.uniq <- subsetByOverlaps(WT, KO, invert = TRUE)
KO.uniq <- subsetByOverlaps(KO, WT, invert = TRUE)

# build objects with the numbers of sites in the subsets
venn <- Venn(SetNames=c("WT", "KO"), 
    Weight=c(
        '10'=length(WT.uniq), 
        '11'=length(subset), 
        '01'=length(KO.uniq)
    )
)

# plot Venn Diagram
pdf("venn_WT_KO.pdf")
plot(venn)
dev.off()
