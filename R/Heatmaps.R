library(DESeq2)
library(pheatmap)

##Data Processing

#Reading in counts, removing X from column names
countdata = read.table("../counts.cts(1)", sep = '\t')
#countdata = read.csv("/Users/elikond/Downloads/counts.cts", sep = '\t')
temp_col_names <- sub('.', '', colnames(countdata)[-c(1)])
new_col_names <- c("Genes",temp_col_names)
colnames(countdata) <- new_col_names

#Removing duplicated and setting row names as genes
countdata <- countdata[!duplicated(countdata$Genes), ]
countdata <- na.omit(countdata)
rownames(countdata) <- countdata$Genes

## Removes Gene column and patients that don't have engraftment data
countdata <- countdata[-c(1,10,12,14,15,17,21,27,28)]

## Reading in metadata, matching patient labels with labels in count data 
coldata = read.csv("../Engraftment - 1.csv", sep = ",")
#coldata = read.csv("/Users/elikond/Downloads/Engraftment - 1.csv", sep = ",")
rownames(coldata) <- sub('-', '.', coldata$Patient)
coldata <- coldata[-c(1)]

#Reading in signature, setting row names as genes
sig9 = read.csv("../sig9.csv", sep = ',')
#sig9 = read.csv("/Users/elikond/Downloads/sig9.csv", sep = ',')
rownames(sig9) <- sig9[,1]
sig9 <- sig9[-c(1)]

## Beginning analysis
# Create DESEQ2 Object
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ Engraftment)

#Keeping genes who's total counts sum to more than one
keep <- rowSums(counts(ddsMat)) > 1
processed_ddsMat <- ddsMat[keep,]

#Apply DESeq and vst (variance stabilizer)
deseq_ob <- DESeq(processed_ddsMat)
vst_ob <- vst(deseq_ob)
res <- results(deseq_ob)
res_df <- as.data.frame(res)

#Gives indeces of rows in deseq_ob with the 20 highest means
select <- order(rowMeans(counts(deseq_ob,normalized=TRUE)), decreasing=TRUE)[1:20]

#Rows = genes, columns = patients 
vst_counts <- assay(vst_ob)
vst_counts_df <- as.data.frame(vst_counts)

#Keeping genes that are in signature
counts_sig9_subset <- vst_counts_df[rownames(vst_counts_df) %in% rownames(sig9), ]

#Making heatmap of genes in sig9)
pheatmap(counts_sig9_subset, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = coldata, fontsize_row = 6.5)

#Making heatmap of genes with highest means
select_counts <- vst_counts_df[select,]
pheatmap(select_counts, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = coldata, fontsize_row = 6.5)

#Making heatmap of genes with highest means that are in sig9
select_sig9_counts <- select_counts[rownames(select_counts) %in% rownames(sig9), ]
pheatmap(select_sig9_counts, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", annotation_col = coldata, fontsize_row = 6.5)

# Consider if we need lfc shrink or not
#lfc_res <- lfcShrink(deseq_ob, coef="Engraftment_l_vs_h", type="apeglm")
#res_df <- as.data.frame(lfc_res)
#res_df

#sum(res$padj < 0.1, na.rm=TRUE)

#Wilcoxon pairwise test on adj pvalue of genes in signature
res_pvals <- res_df[,c('padj','pvalue')]
res_data_subset <- res_pvals[rownames(res_pvals) %in% rownames(sig9), ]
pairwise.wilcox.test(res_data_subset$padj, rownames(res_data_subset),p.adjust.method = "BH", exact=FALSE)
