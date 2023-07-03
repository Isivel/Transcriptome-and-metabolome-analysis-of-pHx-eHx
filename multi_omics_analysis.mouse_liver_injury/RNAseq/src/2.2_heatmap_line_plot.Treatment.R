
##########################################################
########
#### heatmap plot for treatment
load('runtime/2.DEG/res_DEGs.condition.rda')
load('runtime/1.data/fpkm_tab.rdata')
load('runtime/1.data/anno_tab.rdata')


contrasts <- rbind(c(0, 1, 0, 0, 0, 1, 0, 0), # eHx vs pHx in 18
                   c(0, 1, 0, 0, 0, 0, 1, 0), # eHx vs pHx in 36
                   c(0, 1, 0, 0, 0, 0, 0, 1) # eHx vs pHx in 72
)
contrasts <- data.frame(contrast = apply(contrasts, 1, function(x){return(paste0(x, collapse = ""))}), stringsAsFactors = F)
contrasts$name <- c('eHx-pHx:18', 'eHx-pHx:36', 'eHx-pHx:72')

res_DEGs.condition$contrast_alias <- contrasts$name[match(res_DEGs.condition$contrast, contrasts$contrast)]

up_tap <- as.data.frame(table(res_DEGs.condition$contrast_alias[res_DEGs.condition$padj < 0.001 & res_DEGs.condition$log2FoldChange > 1]))
contrasts$up <- up_tap$Freq[match(contrasts$name, up_tap$Var1)]

down_tap <- as.data.frame(table(res_DEGs.condition$contrast_alias[res_DEGs.condition$padj < 0.001 & res_DEGs.condition$log2FoldChange < -1]))
contrasts$down <- down_tap$Freq[match(contrasts$name, down_tap$Var1)]

contrasts$sum <- contrasts$up + contrasts$down
write.table(contrasts, col.names = T, row.names = F, quote = F, sep = "\t", file = 'runtime/2.DEG/DEG_stat.condition.tsv')

res_DEGs.condition$log2FoldChange <- signif(res_DEGs.condition$log2FoldChange, digits = 3)
res_DEGs.condition$stat <- signif(res_DEGs.condition$stat, digits = 3)
res_DEGs.condition$lfcSE <- signif(res_DEGs.condition$lfcSE, digits = 3)
res_DEGs.condition$pvalue <- signif(res_DEGs.condition$pvalue, digits = 3)
res_DEGs.condition$padj <- signif(res_DEGs.condition$padj, digits = 3)
write.table(res_DEGs.condition, col.names = T, row.names = F, quote = F, sep = "\t", file = 'runtime/2.DEG/DEG.condition.tsv')





fpkm_tab <- log2(fpkm_tab + 1)

## number of DEGs is 3729
load('runtime/2.DEG/res_DEGs.condition.rda')
gene.DE <- unique(res_DEGs.condition$geneID[which(res_DEGs.condition$padj < 0.001 &
                                                    abs(res_DEGs.condition$log2FoldChange) > 1)])

annotation_col = data.frame(ID = colnames(fpkm_tab))
annotation_col$Group <- gsub('(^[^-]+)-.*$', '\\1', annotation_col$ID, perl = T)

annotation_col$Time <- gsub('^[^-]+-(\\d+)-.*$', '\\1', annotation_col$ID, perl = T)
annotation_col$Time[grep('sham', annotation_col$Time)] <- 0

rownames(annotation_col) <- annotation_col$ID
annotation_col$ID <- NULL

ann_colors = list(
  Group = c("sham"="steelblue","eHx"="red", "pHx"="#FF7F00"),
  Time = c('0' = 'navy', '18' = '#ffcccb', '36' = 'salmon', '72' = 'darkred')
)
sample_order<-order(annotation_col$Group, annotation_col$Time)
pdf("runtime/2.DEG/heatmap.treatment.pdf", width = 7.5, height = 7)
out_treatment <- pheatmap::pheatmap(fpkm_tab[gene.DE,],
                                    clustering_method = "ward.D2",
                                    clustering_distance_rows = 'correlation',
                                    show_rownames =F,show_colnames =F,fontsize_row=7,
                                    annotation_col = annotation_col,cluster_cols = F,scale = 'row',
                                    annotation_colors = ann_colors,color = colorRampPalette(c("blue", "white", "firebrick3"))(100)) 
dev.off()

save(out_treatment, file = "runtime/2.DEG/pheatmap_treatment.Rdata")

cluster_tab.treatment <- sort(cutree(out_treatment$tree_row, k = 5))
annotation_row <- data.frame(ID = rownames(fpkm_tab), stringsAsFactors = F)
annotation_row$Cluster <- cluster_tab.treatment[match(annotation_row$ID, names(cluster_tab.treatment))]
annotation_row$Cluster <- paste0('Cluster', annotation_row$Cluster)
rownames(annotation_row) <- annotation_row$ID
annotation_row <- annotation_row[annotation_row$Cluster != 'ClusterNA',]
tmp_tab <- annotation_row
tmp_tab$geneSymbol <- res_DEGs.condition$geneSymbol[match(tmp_tab$ID, rownames(res_DEGs.condition))]
write.table(tmp_tab, col.names = T, row.names = F, sep = "\t", quote = F,
            file = "runtime/2.DEG/heatmap.treatment.cluster.txt")
annotation_row$ID <- NULL

require(randomcoloR)
n <- 5
palette <- distinctColorPalette(n)
ann_colors = list(
  Group = c("sham"="steelblue","eHx"="red", "pHx"="#FF7F00"),
  Time = c('0' = 'navy', '18' = '#ffcccb', '36' = 'salmon', '72' = 'darkred'),
  Cluster = c('Cluster1' = palette[1], 'Cluster2' = palette[2], 'Cluster3' = palette[3], 'Cluster4' = palette[4],
              'Cluster5' = palette[5])
)

pdf("runtime/2.DEG/heatmap.treatment.withRowAnno.pdf", width = 7.5, height = 7)
out_treatment.new <- pheatmap::pheatmap(fpkm_tab[gene.DE,],
                                        clustering_method = "ward.D2",
                                        clustering_distance_rows = 'correlation',
                                        show_rownames =F,show_colnames =F,fontsize_row=7,
                                        annotation_col = annotation_col,
                                        annotation_row = annotation_row,
                                        cluster_cols = F,scale = 'row',
                                        annotation_colors = ann_colors,color = colorRampPalette(c("blue", "white", "firebrick3"))(100)) 
dev.off()


## trend graph
dat_expr <- fpkm_tab[names(cluster_tab.treatment),]
dat_expr <- cbind(dat_expr[,1:3], dat_expr)
colnames(dat_expr)[4:6] <- paste0(colnames(dat_expr)[4:6], "_2")
dat_expr <- as.data.frame(t(dat_expr))
dat_expr$ID <- rownames(dat_expr)
dat_tmp <- melt(dat_expr)
dat_tmp$condition <- gsub('(^[^-]+)-.*$', '\\1', dat_tmp$ID, perl = T)
dat_tmp$condition[grepl('sham', dat_tmp$ID)] <- 'pHx'
dat_tmp$condition[grepl('_2$', dat_tmp$ID)] <- 'eHx'
dat_tmp$Time <- gsub('^[^-]+-(\\d+)-.*$', '\\1', dat_tmp$ID, perl = T)
dat_tmp$Time[grep('sham', dat_tmp$Time)] <- 0

df_plot <- aggregate(value ~ Time + condition + variable, data = dat_tmp, mean)

df_plot$Time <- as.numeric(df_plot$Time)
df_plot$cluster <- cluster_tab.treatment[match(df_plot$variable, names(cluster_tab.treatment))]
df_plot$cluster <- paste0('Cluster ', df_plot$cluster)

df_plot$ID <- paste0(df_plot$condition, '_', df_plot$variable)

plot_list <- list()
for(clu in unique(df_plot$cluster)){
  df_plot.sub <- df_plot[df_plot$cluster == clu, ]
  #df_plot.sub$ID <- factor(df_plot.sub$ID)
  
  line_plot <- aggregate(value ~ Time + condition, data = df_plot.sub, mean)
  line_plot$cluster <- unique(df_plot.sub$cluster)
  line_plot$variable <- paste0(line_plot$condition, 'Col')
  line_plot$ID <- line_plot$variable
  
  df_plot.sub <- rbind(df_plot.sub, line_plot[, c('Time', 'condition','variable', 'value', 'cluster',  'ID')])
  
  #df_plot.sub <- df_plot.sub[grep('Col', df_plot.sub$ID),]  
  df_plot.sub$ID <- factor(df_plot.sub$ID)
  
  cols <- rep('grey', nlevels(df_plot.sub$ID))
  cols[grep('pHxCol',levels(df_plot.sub$ID), perl = T)] <- 'steelblue'
  cols[grep('eHxCol',levels(df_plot.sub$ID), perl = T)] <- 'salmon'
  
  df_plot.sub$al <- 0.6
  df_plot.sub$al[grep('Col$', df_plot.sub$ID, perl = T)] <- 1
  
  
  
  plot_list[[clu]] <- ggplot(data = df_plot.sub, aes(x = Time, y = value, color = ID)) +
    geom_line(aes(alpha = al)) +
    scale_color_manual(values = cols) +
    #scale_alpha_discrete(range = c(0.6, 1)) +
    scale_x_continuous(breaks = c(0, 18, 36, 72), labels = c('0','18','36','72')) +
    labs(y = 'log2(FPKM + 1)', title = clu) +
    theme(
      plot.title = element_text(hjust = .5, vjust = .5, size = 18),
      axis.text.x = element_text(size = 13),
      axis.title = element_text(size = 15),
      #axis.title.x = element_blank(),
      #legend.title = element_blank(),
      legend.position = 'none',
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line()
    )
}

require(ggpubr)
pdf("runtime/2.DEG/DEG_expr.line.treatment.pdf",width = 17,height = 9, onefile = T)
ggarrange(plotlist = plot_list,ncol = 3, nrow = 2)
dev.off()

save(df_plot, file = 'runtime/2.DEG/data_plot_line.treatment.Rdata')

## save data frame with cluster information
load('runtime/1.data/anno_tab.rdata')
res_cluster <- data.frame(geneID = rownames(fpkm_tab), stringsAsFactors = F)
res_cluster$Cluster <- paste0('Cluster', cluster_tab.treatment[match(res_cluster$geneID, names(cluster_tab.treatment))])
res_cluster$geneSymbol <- anno_tab$Name[match(res_cluster$geneID, anno_tab$Gene.ID)]
write.table(res_cluster, col.names = T, row.names = F, sep = "\t", quote = F, 
            file = "runtime/2.DEG/GEgene_cluster.treatment.tsv")
