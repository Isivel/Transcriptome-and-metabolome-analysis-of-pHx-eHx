
library(openxlsx)
library(genefilter)
library(DESeq2)
library(ggplot2)
#library(devtools)
#install_github("wjawaid/enrichR")
library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(stringr)
library(reshape2)

outDir <- './runtime/1.data'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

expr_tab <- read.xlsx('Expression.xlsx', check.names = F)
colnames(expr_tab) <- gsub('^pH-', 'pHx-', colnames(expr_tab), perl = T)

anno_tab <- read.xlsx('Annotation.xlsx', check.names = F)
save(anno_tab, file = 'runtime/1.data/anno_tab.rdata')

ind_count <- c(1, grep('read.count', colnames(expr_tab)))
count_tab <- expr_tab[, ind_count]
colnames(count_tab) <- gsub('\\:read.count', '', colnames(count_tab))
rownames(count_tab) <- count_tab$Gene_ID
count_tab$Gene_ID <- NULL
save(count_tab, file = 'runtime/1.data/count_tab.rdata')

ind_fpkm <- c(1, grep('fpkm', colnames(expr_tab)))
fpkm_tab <- expr_tab[, ind_fpkm]
colnames(fpkm_tab) <- gsub('\\:fpkm', '', colnames(fpkm_tab))
rownames(fpkm_tab) <- fpkm_tab$Gene_ID
fpkm_tab$Gene_ID <- NULL

save(fpkm_tab, file = 'runtime/1.data/fpkm_tab.rdata')

#### sample cluster plot
fpkm_tab.sub <- varFilter(as.matrix(fpkm_tab), var.func=IQR,
                          var.cutoff=0.5, filterByQuantile=TRUE)

cor_mat <- cor(fpkm_tab.sub, method = 'pearson')
clusters <- hclust(dist(1 - cor_mat))
pdf(file = 'mRNA_hclust.pdf', width = 6, height = 6)
plot(clusters, hang = -1, xlab = NA, sub = NA)
dev.off()



