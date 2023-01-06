library(DESeq2)
library(tximeta)
library(readr)
library(tximeta)
library(apeglm)
library(kohonen)
library(tibble)
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(grid)
library(plotly)

countdata = read_tsv("/Users/elikond/Downloads/counts.cts", col_names = TRUE)
countdata <- as.data.frame(countdata)
countdata <- countdata[!duplicated(countdata$Genes), ]
countdata <- na.omit(countdata)
rownames(countdata) <- countdata$Genes
countdata <- countdata[-c(1,10,12,14,15,17,21,27,28)]

coldata = read_tsv("/Users/elikond/Downloads/Engraftment_Sheet2.tsv", col_names = TRUE)

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ Engraftment)

keep <- rowSums(counts(ddsMat)) > 1
ddsMat <- ddsMat[keep,]

ddsMat <- DESeq(ddsMat)
res <- results(ddsMat)

resultsNames(ddsMat)
res <- lfcShrink(ddsMat, coef="Engraftment_l_vs_h", type="apeglm")

df_res <- as.data.frame(res)
df_res <- df_res[-c(1:3)]
setnames(df_res, old = c('pvalue', 'padj'), new = c('pvalue_bulk', 'padj_bulk'))

my_sig9 = read_csv("/Users/elikond/Downloads/sig9.csv", col_names = TRUE)
my_sig9 <- as.data.frame(my_sig9)
rownames(my_sig9) <- my_sig9[,1]
my_sig9 <- my_sig9[-c(1,3:5)]
setnames(my_sig9, old = c('p_val', 'p_val_adj'), new = c('pvalue_sig', 'padj_sig'))

df2 <- merge(my_sig9, df_res, by=0)
rownames(df2) <- df2[,1]
df2 <- df2[-c(1)]

data <- apply(df2, MARGIN = 1, function(x){x/mean(x)})
data <- t(data)

data <- as.data.frame(data)

new_df <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))
for (x in colnames(data)){
  new_df[,x] = qnorm(data[[x]])
}


new_df <- new_df[-c(1:4,8)]
setnames(new_df, old = c('pvalue_sig', 'padj_sig', 'pvalue_bulk'), new = c('zscore_sig','zscore_adj_sig','zscore_bulk'))
rownames(new_df) <- rownames(data)

plot_matrix <- as.matrix(new_df)
p <- plot_ly(x=colnames(plot_matrix), y=rownames(plot_matrix), z = plot_matrix, type = "heatmap") %>%
  layout(margin = list(l=120))
print(p)
