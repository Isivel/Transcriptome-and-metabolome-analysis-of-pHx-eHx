
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

outDir <- './runtime/2.DEG'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

#### call DEGs
load('runtime/1.data/count_tab.rdata')
load('runtime/1.data/fpkm_tab.rdata')

## remove low expression feature
group <- gsub('-\\d$', '', colnames(count_tab), perl = T)

keep <- apply(count_tab, 1, function(x){
  x <- as.numeric(x)
  df_tmp <- data.frame(group, expr = x, stringsAsFactors = F)
  mean_g <- aggregate(expr~group, data = df_tmp, mean)
  any(mean_g$expr >= 3)
})

count_tab <- count_tab[keep,]
#dim(count_tab)
#[1] 14505    21
write.table(data.frame(gene = rownames(count_tab), count_tab), row.names = F,
            col.names = T, sep = "\t", quote = F, 
            file = 'count_table.filtered.txt')

## add repeat control sample for avoiding 'Model matrix not full rank' in DESeq2
count_tab <- cbind(count_tab[, 1:3], count_tab)
colnames(count_tab)[4:6] <- paste0(colnames(count_tab)[4:6], "_2")

## make meta data
col_tab <- data.frame(ID = colnames(count_tab),
                      stringsAsFactors = F)
col_tab$condition <- gsub('(^[^-]+)-.*$', '\\1', col_tab$ID, perl = T)
col_tab$condition[grepl('sham', col_tab$ID)] <- 'pHx'
col_tab$condition[grepl('_2$', col_tab$ID)] <- 'eHx'
col_tab$condition <- factor(col_tab$condition, levels = c('pHx', 'eHx'))

col_tab$time <- gsub('^[^-]+-(\\d+)-.*$', '\\1', col_tab$ID, perl = T)
col_tab$time[grep('sham', col_tab$time)] <- 0
col_tab$time <- factor(col_tab$time)

dd <- DESeqDataSetFromMatrix(
  countData = count_tab,
  colData = col_tab,
  design = ~ condition + time + condition:time
)

#design(dd) = as.formula(paste0('~ ', variable_to_design))
dds <- DESeq(dd)

save(dds, file = 'runtime/2.DEG/dds.rdata')


##### extract DEG result in time point
contrasts <- rbind(c(0, 0, 1, 0, 0, 0, 0, 0), # 18 vs 0 in pHx
                   c(0, 0, 0, 1, 0, 0, 0, 0), # 36 vs 0 in pHx
                   c(0, 0, 0, 0, 1, 0, 0, 0), # 72 vs 0 in pHx
                   c(0, 0, 1, 0, 0, 1, 0, 0), # 18 vs 0 in eHx
                   c(0, 0, 0, 1, 0, 0, 1, 0), # 36 vs 0 in eHx
                   c(0, 0, 0, 0, 1, 0, 0, 1), # 72 vs 0 in eHx
                   c(0, 0, -1, 1, 0, 0, 0, 0), # 36 vs 18 in pHx
                   c(0, 0, -1, 0, 1, 0, 0, 0), # 72 vs 18 in pHx
                   c(0, 0, -1, 1, 0, -1, 1, 0), # 36 vs 18 in eHx
                   c(0, 0, -1, 0, 1, -1, 0, 1), # 72 vs 18 in eHx
                   c(0, 0, 0, -1, 1, 0, 0, 0), # 72 vs 36 in pHx
                   c(0, 0, 0, -1, 1, 0, -1, 1) # 72 vs 36 in eHx
                   )

res_DEGs.time <- NULL
for(i in 1:nrow(contrasts)){
  contrast <- contrasts[i,]
  res <- as.data.frame(results(dds, contrast = contrast))
  res$contrast <- paste0(contrast, collapse = "")
  res$geneID <- rownames(res)
  res_DEGs.time <- rbind(res_DEGs.time, res)
}

load('runtime/1.data/anno_tab.rdata')
res_DEGs.time$geneSymbol <- anno_tab$Name[match(res_DEGs.time$geneID, anno_tab$Gene.ID)]
save(res_DEGs.time, file = 'runtime/2.DEG/res_DEGs.time.rda')

res_DEGs.time$baseMean <- round(res_DEGs.time$baseMean, digits = 4)
res_DEGs.time$log2FoldChange <- signif(res_DEGs.time$log2FoldChange, digits = 4)
res_DEGs.time$lfcSE <- signif(res_DEGs.time$lfcSE, digits = 4)
res_DEGs.time$stat <- signif(res_DEGs.time$stat, digits = 4)
res_DEGs.time$pvalue <- signif(res_DEGs.time$pvalue, digits = 4)
res_DEGs.time$padj <- signif(res_DEGs.time$padj, digits = 4)
write.table(res_DEGs.time, col.names = T, row.names = F, quote = F, sep = "\t",
            file = 'runtime/2.DEG/res_DEGs.time.tsv')


contrasts <- rbind(c(0, 1, 0, 0, 0, 1, 0, 0), # eHx vs pHx in 18
                   c(0, 1, 0, 0, 0, 0, 1, 0), # eHx vs pHx in 36
                   c(0, 1, 0, 0, 0, 0, 0, 1) # eHx vs pHx in 72
)

res_DEGs.condition <- NULL
for(i in 1:nrow(contrasts)){
  contrast <- contrasts[i,]
  res <- as.data.frame(results(dds, contrast = contrast))
  res$contrast <- paste0(contrast, collapse = "")
  res$geneID <- rownames(res)
  res_DEGs.condition <- rbind(res_DEGs.condition, res)
}

res_DEGs.condition$geneSymbol <- anno_tab$Name[match(res_DEGs.condition$geneID, anno_tab$Gene.ID)]
save(res_DEGs.condition, file = 'runtime/2.DEG/res_DEGs.condition.rda')

res_DEGs.condition$baseMean <- round(res_DEGs.condition$baseMean, digits = 4)
res_DEGs.condition$log2FoldChange <- signif(res_DEGs.condition$log2FoldChange, digits = 4)
res_DEGs.condition$lfcSE <- signif(res_DEGs.condition$lfcSE, digits = 4)
res_DEGs.condition$stat <- signif(res_DEGs.condition$stat, digits = 4)
res_DEGs.condition$pvalue <- signif(res_DEGs.condition$pvalue, digits = 4)
res_DEGs.condition$padj <- signif(res_DEGs.condition$padj, digits = 4)
write.table(res_DEGs.condition, col.names = T, row.names = F, quote = F, sep = "\t",
            file = 'runtime/2.DEG/res_DEGs.condition.tsv')











