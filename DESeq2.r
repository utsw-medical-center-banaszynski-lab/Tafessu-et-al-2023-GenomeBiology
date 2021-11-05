library(DESeq2)

expr_readcount_combined <- read.table("combined_counts_WT_KO.txt", header = T, row.names= 1, sep="\t")

# create sample information table
col_data <- data.frame(library = colnames(expr_readcount_combined),
                       sample = c(rep("KO",2),
                                  rep("WT",2)))

col_data$sample <- factor(col_data$sample, 
                          levels = c("KO", "WT"))

# create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = expr_readcount_combined,
                              colData = col_data,
                              design = ~ sample)

# perform DE detection
dds <- DESeq(dds)
res <- results(dds, 
               contrast = c("sample", "KO", "WT"),
               alpha = 0.05)

# output results
write.table(cbind(gene_name = rownames(as.data.frame(res)), 
                  as.data.frame(res)), 
            file = "DESeq2_WT_KO.txt",
            sep = "\t",
            quote = F,
            row.names = F,
            col.names =T)