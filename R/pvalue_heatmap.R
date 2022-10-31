install.packages("apeglm")

BiocManager::install("apeglm")

library(DESeq2)
library(tximeta)
library(apeglm)
library(kohonen)
library(tibble)
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(grid)
library(plotly)
library("pheatmap")


countdata = read.table("../counts(1).cts", sep = '\t')


colnames(countdata) <- countdata[1,]
countdata <- countdata[-c(1), ]




countdata <- countdata[!duplicated(countdata$Genes), ]
countdata <- na.omit(countdata)
rownames(countdata) <- countdata$Genes
## Removes patients we didnt scRNAseq
countdata <- countdata[-c(1,10,12,14,15,17,21,27,28)]

## convert df to numeric
countdata <- mutate_all(countdata, function(x) as.numeric(as.character(x))) 

## load in metadata
coldata = read.csv("Engraftment_Sheet2.tsv", sep = "\t")


## create deseq2 object
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ Engraftment)



keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]

ddsMat <- DESeq(ddsMat)
vst_ob <- vst(ddsMat)
res <- results(ddsMat)

resultsNames(ddsMat)


vst_ob

#### HEATMAP WES CODE STARTS HERE



select <- order(rowMeans(counts(ddsMat,normalized=TRUE)), decreasing=TRUE)[1:20]

df8 <- assay(vst_ob)

df_3 <- colData(ddsMat)[,"Engraftment"]

supData <- rownames(sig9)

supData <- supData[-c(16)]

df8 <- as.data.frame(df8)

pheatmap(df8[c(paste(supData)),], cluster_rows = FALSE, cluster_cols = TRUE, scale = "row")


#### HEATMAP WES CODE ENDS HERE


# Consider if we need lfc shrink or not
lfc_res <- lfcShrink(ddsMat, coef="Engraftment_l_vs_h", type="apeglm")
res_df <- as.data.frame(lfc_res)
res_df

# removing cols for E
res_df <- res_df[-c(1:3,5)]

sig9 = read_csv("../sig9.csv", col_names = TRUE)
sig9 <- as.data.frame(sig9)

# renaming rownames
rownames(sig9) <- sig9[,1]

sig9

# removing extra Calcs
sig9 <- sig9[-c(1,3:6)]
setnames(sig9, old = c('p_val'), new = c('pvalue_sig'))

df1 <- merge(sig9, res_df, by=0, all=TRUE)
df1 <- na.omit(df1)
rownames(df1) <- df1[,1]
df1 <- df1[-c(1,2)]

final_data <- as.matrix(df1)
#data <- apply(data, MARGIN = 1, function(x){x/mean(x)})

p <- plot_ly(x=colnames(final_data), y=rownames(final_data), z = final_data, type = "heatmap") %>%
  layout(margin = list(l=120))
print(p)

zscore_df <- data.frame(matrix(ncol = ncol(df1), nrow = nrow(df1)))
zscore_df$zscore = qnorm(df1$pvalue)

zscore_df <- zscore_df[-c(1)]
rownames(zscore_df) <- rownames(df1)
pairwise.wilcox.test(zscore_df$zscore, rownames(zscore_df),p.adjust.method = "BH")
