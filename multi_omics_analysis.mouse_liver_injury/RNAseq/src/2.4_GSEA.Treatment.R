
library(ggplot2)
library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(stringr)
library(reshape2)

##### GSEA
dotplot_ylab <- function(x,showCategory = 20, wi = 50, term_cex = 8, term = NA){
  x.df <- as.data.frame(x)
  
  if(nrow(x.df)>showCategory){
    x.df <- x.df[1:showCategory,]
  }
  x.df$GR <- as.numeric(gsub("^(\\d+)\\/\\d+$", "\\1",x.df$GeneRatio)) / as.numeric(gsub("^\\d+\\/(\\d+)$", "\\1",x.df$GeneRatio))
  #x.df$oddRatios <-  unlist(lapply(as.list(x.df$oddRatios), function(x) eval(parse(text=x))))
  x.df$Description <- factor(x.df$Description,levels = x.df$Description[order(x.df$GR,decreasing = F)])
  
  #x.df <- x.df[order(x.df$Count,decreasing = F),]
  p<-ggplot(x.df, aes(x=GR, y=Description,color=pvalue)) +
    geom_point(aes(size = Count))+scale_color_gradient(low="red", high="blue") +
    labs(x = "GeneRatio", title = term) +
    scale_y_discrete(labels=function(x) str_wrap(x,width=wi))+
    theme(axis.text.y = element_text( size = term_cex),
          plot.title = element_text(hjust = .5)) +
    guides(color = guide_legend(order = 0),
           size = guide_legend(order = 1))
}
### get gene list in interest module
## cluster high-low-mid-high

cluster_tab <- read.delim('runtime/2.DEG/GEgene_cluster.treatment.tsv', stringsAsFactors = F, check.names = F)

#interest cluster
interest_cluster <- c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5')
for(clu in interest_cluster){
  modProbes = unique(as.character(cluster_tab$geneID[which(cluster_tab$Cluster %in% clu)]))
  
  #### enrichment analysis
  
  go_BP <- enrichGO(gene = modProbes, OrgDb = 'org.Mm.eg.db',
                    keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.2)
  go_MF <- enrichGO(gene = modProbes, OrgDb = 'org.Mm.eg.db',
                    keyType = 'ENSEMBL', ont = 'MF', pvalueCutoff = 0.2)
  go_CC <- enrichGO(gene = modProbes, OrgDb = 'org.Mm.eg.db',
                    keyType = 'ENSEMBL', ont = 'CC', pvalueCutoff = 0.2)
  
  dot_BP <- dotplot_ylab(go_BP, showCategory=20, term = 'GO_BP')
  dot_MF <- dotplot_ylab(go_MF, showCategory=20, term = 'GO_MF')
  dot_CC <- dotplot_ylab(go_CC, showCategory=20, term = 'GO_CC')
  
  plotDir <- paste0("runtime/2.DEG/GSEA/Treatment/", clu)
  if(!dir.exists(plotDir)){
    dir.create(plotDir, recursive = TRUE)
  }
  
  ggsave(filename = paste0(plotDir, '/GO_BP','.pdf'), 
         width = 8, height = 8,plot = dot_BP)
  
  ggsave(filename = paste0(plotDir, '/GO_CC','.pdf'), 
         width = 8, height = 8,plot = dot_CC)
  
  ggsave(filename = paste0(plotDir, '/GO_MF','.pdf'),  
         width = 8, height = 8,plot = dot_MF)
  
  
  write.table(go_BP, col.names = T, row.names = F,
              sep = "\t", quote = F, 
              file = paste0(plotDir,'/GO_BP.enrichment.txt'))
  write.table(go_CC, col.names = T, row.names = F,
              sep = "\t", quote = F, 
              file = paste0(plotDir,'/GO_CC.enrichment.txt'))
  write.table(go_MF, col.names = T, row.names = F,
              sep = "\t", quote = F, 
              file = paste0(plotDir,'/GO_MF.enrichment.txt'))
  
  
  modProbes_entrez <- bitr(modProbes, fromType = 'ENSEMBL',
                           toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db')
  kk <- enrichKEGG(gene=modProbes_entrez$ENTREZID, organism='mmu', pvalueCutoff = 1,
                   qvalueCutoff = 1)
  
  dot_KK <- dotplot_ylab(kk, showCategory=20, term = 'KEGG')
  
  ggsave(filename = paste0(plotDir, '/KEGG','.pdf'),  
         width = 8, height = 8,plot = dot_KK)
  write.table(kk, col.names = T, row.names = F,
              sep = "\t", quote = F, 
              file = paste0(plotDir,'/KEGG.enrichment.txt'))
}




