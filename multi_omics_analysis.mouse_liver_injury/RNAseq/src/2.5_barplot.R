
library(ggplot2)
library(reshape2)


load('runtime/1.data/fpkm_tab.rdata')
load('runtime/1.data/anno_tab.rdata')
fpkm_tab <- log2(fpkm_tab + 1)

cluster_tab <- read.delim('runtime/2.DEG/GEgene_cluster.time.tsv', stringsAsFactors = F, check.names = F)

load('runtime/2.DEG/res_DEGs.time.rda')

plotDir <- 'runtime/2.DEG/boxplot/Time/'
if(!dir.exists(plotDir)){
  dir.create(plotDir, recursive = T)
}

clusters <- c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4',
             'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8')

for(clu in clusters){
  #clu <- 'Cluster1'
  res_DEGs.time.sub <- res_DEGs.time[res_DEGs.time$geneID %in% cluster_tab$geneID[which(cluster_tab$Cluster == clu)],]
  res_DEGs.time.sub <- res_DEGs.time.sub[order(res_DEGs.time.sub$pvalue),]
  
  genes <- unique(res_DEGs.time.sub$geneID)[1:10]
  
  expr_tab.gene <- as.data.frame(t(fpkm_tab[genes, ]))
  expr_tab.gene$sample <- rownames(expr_tab.gene)
  
  dat_tmp <- melt(expr_tab.gene)
  
  dat_tmp$treatment <- gsub('(^[^-]+)-.*$', '\\1', dat_tmp$sample, perl = T)
  
  dat_tmp$Time <- gsub('^[^-]+-(\\d+)-.*$', '\\1', dat_tmp$sample, perl = T)
  dat_tmp$Time[grep('sham', dat_tmp$Time)] <- 0
  
  dat_tmp$gene <- anno_tab$Name[match(dat_tmp$variable, anno_tab$Gene.ID)]
  dat_tmp$gene <- factor(dat_tmp$gene, levels = anno_tab$Name[match(genes, anno_tab$Gene.ID)])
  
  boxplot <- ggplot(data = dat_tmp, aes(x = Time, y = value)) +
    geom_boxplot(aes(fill = treatment)) +
    labs(y = "log2(FPKM + 1)", fill = 'Group') +
    scale_fill_manual(values = c('sham'='grey', 'pHx' = 'steelblue', 'eHx' = 'salmon')) +
    theme_bw() +
    theme(
      strip.text.x = element_text(size = 12),
      legend.title = element_blank()
    ) + facet_wrap(~gene, ncol = 5, scales = 'free')
  
  ggsave(plot = boxplot, width = 17, height = 7, filename = paste0('runtime/2.DEG/boxplot/Time/', clu, '.pdf'))
  
  
}


cluster_tab <- read.delim('runtime/2.DEG/GEgene_cluster.treatment.tsv', stringsAsFactors = F, check.names = F)

load('runtime/2.DEG/res_DEGs.condition.rda')

plotDir <- 'runtime/2.DEG/boxplot/Treatment/'
if(!dir.exists(plotDir)){
  dir.create(plotDir, recursive = T)
}

clusters <- c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4',
              'Cluster5')

for(clu in clusters){
  #clu <- 'Cluster1'
  res_DEGs.condition.sub <- res_DEGs.condition[res_DEGs.condition$geneID %in% cluster_tab$geneID[which(cluster_tab$Cluster == clu)],]
  res_DEGs.condition.sub <- res_DEGs.condition.sub[order(res_DEGs.condition.sub$pvalue),]
  
  genes <- unique(res_DEGs.condition.sub$geneID)[1:10]
  
  expr_tab.gene <- as.data.frame(t(fpkm_tab[genes, ]))
  expr_tab.gene$sample <- rownames(expr_tab.gene)
  
  dat_tmp <- melt(expr_tab.gene)
  
  dat_tmp$treatment <- gsub('(^[^-]+)-.*$', '\\1', dat_tmp$sample, perl = T)
  
  dat_tmp$Time <- gsub('^[^-]+-(\\d+)-.*$', '\\1', dat_tmp$sample, perl = T)
  dat_tmp$Time[grep('sham', dat_tmp$Time)] <- 0
  
  dat_tmp$gene <- anno_tab$Name[match(dat_tmp$variable, anno_tab$Gene.ID)]
  dat_tmp$gene <- factor(dat_tmp$gene, levels = anno_tab$Name[match(genes, anno_tab$Gene.ID)])
  
  boxplot <- ggplot(data = dat_tmp, aes(x = Time, y = value)) +
    geom_boxplot(aes(fill = treatment)) +
    labs(y = "log2(FPKM + 1)", fill = 'Group') +
    scale_fill_manual(values = c('sham'='grey', 'pHx' = 'steelblue', 'eHx' = 'salmon')) +
    theme_bw() +
    theme(
      strip.text.x = element_text(size = 12),
      legend.title = element_blank()
    ) + facet_wrap(~gene, ncol = 5, scales = 'free')
  
  ggsave(plot = boxplot, width = 17, height = 7, filename = paste0('runtime/2.DEG/boxplot/Treatment/', clu, '.pdf'))
  
  
}
