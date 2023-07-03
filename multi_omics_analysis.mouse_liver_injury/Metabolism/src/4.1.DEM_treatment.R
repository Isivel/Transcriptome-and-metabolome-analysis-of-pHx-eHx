
library(ggplot2)
library(reshape2)

outDir <- './runtime/4.DEM.new/Treatment'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

dat.met <- readRDS('runtime/1.data/all_sample.onlyExpr.RDS')
dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')


dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')
dat.met.norm.twA <- cbind(dat.met.norm[, 1:5], dat.met.norm)
colnames(dat.met.norm.twA)[1:5] <- paste0(colnames(dat.met.norm.twA)[1:5], '_2')
df_calDiff <- data.frame(sample = colnames(dat.met.norm.twA), stringsAsFactors = F)
df_calDiff$Group <- gsub('^([^-]+)\\-.*$', '\\1', df_calDiff$sample, perl = T)
df_calDiff$Time <- gsub('^.*\\-(\\d+)\\-.*$', '\\1', df_calDiff$sample, perl = T)
df_calDiff$Time[grep('sham', df_calDiff$Time)] <- '0'
df_calDiff$Group[1:5] <- 'pH'
df_calDiff$Group[6:10] <- 'eHx'


load('runtime/4.DEM.new/res_DEM.rda')

compounds <- unique(res_DEM$Index[which(grepl('eHx-pH', res_DEM$contrast) & res_DEM$q < 0.05)])

plotDir <- 'runtime/4.DEM.new/Treatment'
## plot heatmap
annotation_col = df_calDiff[!grepl('sham.*_2',df_calDiff$sample, perl = T), ]
annotation_col$Group[grepl('sham', annotation_col$sample)] <- 'sham'
rownames(annotation_col) <- annotation_col$sample
annotation_col$sample <- NULL
ann_colors = list(
  Group = c("sham"="steelblue","eHx"="red", "pH"="#FF7F00"),
  Time = c('0' = 'navy', '18' = '#ffcccb', '36' = 'salmon', '72' = 'darkred')
)
sample_order<-order(annotation_col$Group, annotation_col$Time)
pdf("runtime/4.DEM.new/Treatment/heatmap.treatment.pdf", width = 7.5, height = 7.5)
out_treatment <- pheatmap::pheatmap(dat.met.norm[compounds,],
                                    clustering_method = "ward.D2",
                                    clustering_distance_rows = 'correlation',
                                    show_rownames =F,show_colnames =F,fontsize_row=7,
                                    annotation_col = annotation_col,cluster_cols = F,scale = 'row',
                               annotation_colors = ann_colors,color = colorRampPalette(c("blue", "white", "firebrick3"))(100)) 
dev.off()

save(out_treatment, file = "runtime/4.DEM.new/Treatment/pheatmap_Treatment.Rdata")


#############################################################
## extract metabolist with cluster
all_data <- readRDS('runtime/1.data/all_sample.withAnno.RDS')



cluster_tab.treatment <- sort(cutree(out_treatment$tree_row, k = 2))
annotation_row <- data.frame(ID = compounds, stringsAsFactors = F)
annotation_row$Cluster <- cluster_tab.treatment[match(annotation_row$ID, names(cluster_tab.treatment))]
annotation_row$Cluster <- paste0('Cluster', annotation_row$Cluster)
rownames(annotation_row) <- annotation_row$ID
annotation_row$ID <- NULL
require(randomcoloR)
n <- 2
palette <- distinctColorPalette(n)
ann_colors = list(
  Group = c("sham"="steelblue","eHx"="red", "pH"="#FF7F00"),
  Time = c('0' = 'navy', '18' = '#ffcccb', '36' = 'salmon', '72' = 'darkred'),
  Cluster = c('Cluster1' = palette[1], 'Cluster2' = palette[2])
)
pdf("runtime/4.DEM.new/Treatment/heatmap.treatment.withCluster.pdf", width = 7.5, height = 7.5)
out_treatment <- pheatmap::pheatmap(dat.met.norm[compounds,],
                                    clustering_method = "ward.D2",
                                    clustering_distance_rows = 'correlation',
                                    show_rownames =F,show_colnames =F,fontsize_row=7,
                                    annotation_col = annotation_col,
                                    annotation_row = annotation_row,
                                    cluster_cols = F,scale = 'row',
                                    annotation_colors = ann_colors,color = colorRampPalette(c("blue", "white", "firebrick3"))(100)) 
dev.off()


meta_tab.treatment <- sort(cutree(out_treatment$tree_row, k = 2))

compound.list <- all_data$Compounds[which(all_data$Index %in% names(meta_tab.treatment)[meta_tab.treatment == 1])]
write.table(compound.list, row.names = F, col.names = F, quote = F, file = 'runtime/4.DEM.new/Treatment/compound_treatment_cluster1.tsv')

compound_cluster <- all_data[which(all_data$Index %in% names(meta_tab.treatment)[meta_tab.treatment == 1]),]
write.table(compound_cluster, row.names = T, col.names = T, quote = F, sep = "\t",
            file = 'runtime/4.DEM.new/Treatment/treatment_Cluster1.tsv')
df_plot.pie <- as.data.frame(table(compound_cluster$Class.I), stringsAsFactors = F)
df_plot.pie$lab <- paste0(df_plot.pie$Var1, ' (', df_plot.pie$Freq, ')')
df_plot.pie$lab <- factor(df_plot.pie$lab, levels = df_plot.pie$lab[order(df_plot.pie$Freq, decreasing = T)])
## pie
pie1 <- ggplot(data = df_plot.pie, aes(x="", y=Freq, fill=lab)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(
    legend.title = element_blank()
  )
ggsave(plot = pie1, width = 9, height = 6, filename = 'runtime/4.DEM.new/Treatment/compound_Cluster1.pie.pdf')
## trend line
df_tmp <- as.data.frame(t(dat.met.norm.twA[names(meta_tab.treatment)[which(meta_tab.treatment == 1)],]))
df_tmp$sample <- rownames(df_tmp)
df_plot.tmp <- melt(df_tmp)
df_plot.tmp$Time <- df_calDiff$Time[match(df_plot.tmp$sample, df_calDiff$sample)]
df_plot.tmp$Group <- df_calDiff$Group[match(df_plot.tmp$sample, df_calDiff$sample)]
df_plot.tmp$ID <- paste0(df_plot.tmp$variable, "_", df_plot.tmp$Group)

df_plot.line <- aggregate(value ~ ID + Time, data = df_plot.tmp, median)
df_plot.line$Group <- gsub('^.*_', '', df_plot.line$ID, perl = T)

line_plot <- aggregate(value ~ Time + Group, data = df_plot.line, mean)
line_plot$ID <- paste0(line_plot$Group, 'Col')

df_plot.line <- rbind(df_plot.line, line_plot[, c('ID', 'Time', 'value', 'Group')])

df_plot.line$ID <- factor(df_plot.line$ID)
cols <- rep('grey', nlevels(df_plot.line$ID))
#cols[grep('sham', df_plot.line$sample)] <- 'grey'
cols[grep('pHCol', levels(df_plot.line$ID))] <- 'steelblue'
cols[grep('eHxCol', levels(df_plot.line$ID))] <- 'salmon'

df_plot.line$al <- 0.6
df_plot.line$al[grepl('Col', df_plot.line$ID)] <- 1

df_plot.line$Time <- as.numeric(df_plot.line$Time)

line_chart <- ggplot(data = df_plot.line, aes(x = Time, y = value, color = ID)) +
  geom_line(aes(alpha = al)) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 18, 36, 72), labels = c('0','18','36','72')) +
  labs(y = 'Normalized scaled intensity') +
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

ggsave(plot = line_chart, width = 6.5, height = 6, filename = 'runtime/4.DEM.new/Treatment/Cluster1_linePlot.pdf')


###### cluster2
compound.list <- all_data$Compounds[which(all_data$Index %in% names(meta_tab.treatment)[meta_tab.treatment == 2])]
write.table(compound.list, row.names = F, col.names = F, quote = F, file = 'runtime/4.DEM.new/Treatment/compound_treatment_cluster2.tsv')

compound_cluster <- all_data[which(all_data$Index %in% names(meta_tab.treatment)[meta_tab.treatment == 2]),]
write.table(compound_cluster, row.names = T, col.names = T, quote = F, sep = "\t",
            file = 'runtime/4.DEM.new/Treatment/treatment_Cluster2.tsv')
df_plot.pie <- as.data.frame(table(compound_cluster$Class.I), stringsAsFactors = F)
df_plot.pie$lab <- paste0(df_plot.pie$Var1, ' (', df_plot.pie$Freq, ')')
df_plot.pie$lab <- factor(df_plot.pie$lab, levels = df_plot.pie$lab[order(df_plot.pie$Freq, decreasing = T)])
## pie
pie1 <- ggplot(data = df_plot.pie, aes(x="", y=Freq, fill=lab)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  theme(
    legend.title = element_blank()
  )
ggsave(plot = pie1, width = 9, height = 6, filename = 'runtime/4.DEM.new/Treatment/compound_Cluster2.pie.pdf')
## trend line
df_tmp <- as.data.frame(t(dat.met.norm.twA[names(meta_tab.treatment)[which(meta_tab.treatment == 2)],]))
df_tmp$sample <- rownames(df_tmp)
df_plot.tmp <- melt(df_tmp)
df_plot.tmp$Time <- df_calDiff$Time[match(df_plot.tmp$sample, df_calDiff$sample)]
df_plot.tmp$Group <- df_calDiff$Group[match(df_plot.tmp$sample, df_calDiff$sample)]
df_plot.tmp$ID <- paste0(df_plot.tmp$variable, "_", df_plot.tmp$Group)

df_plot.line <- aggregate(value ~ ID + Time, data = df_plot.tmp, median)
df_plot.line$Group <- gsub('^.*_', '', df_plot.line$ID, perl = T)

line_plot <- aggregate(value ~ Time + Group, data = df_plot.line, mean)
line_plot$ID <- paste0(line_plot$Group, 'Col')

df_plot.line <- rbind(df_plot.line, line_plot[, c('ID', 'Time', 'value', 'Group')])

df_plot.line$ID <- factor(df_plot.line$ID)
cols <- rep('grey', nlevels(df_plot.line$ID))
#cols[grep('sham', df_plot.line$sample)] <- 'grey'
cols[grep('pHCol', levels(df_plot.line$ID))] <- 'steelblue'
cols[grep('eHxCol', levels(df_plot.line$ID))] <- 'salmon'

df_plot.line$al <- 0.6
df_plot.line$al[grepl('Col', df_plot.line$ID)] <- 1

df_plot.line$Time <- as.numeric(df_plot.line$Time)

line_chart <- ggplot(data = df_plot.line, aes(x = Time, y = value, color = ID)) +
  geom_line(aes(alpha = al)) +
  scale_color_manual(values = cols) +
  scale_x_continuous(breaks = c(0, 18, 36, 72), labels = c('0','18','36','72')) +
  labs(y = 'Normalized scaled intensity') +
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

ggsave(plot = line_chart, width = 6.5, height = 6, filename = 'runtime/4.DEM.new/Treatment/Cluster2_linePlot.pdf')





#### boxplot
dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')


### for Cluster1
for(clu in c('1', '2')){
  dat.compound_DE <- all_data[which(all_data$Index %in% names(meta_tab.treatment)[meta_tab.treatment == clu]),]
  res_DEM.sub <- res_DEM[res_DEM$Index %in% dat.compound_DE$Index & grepl('eHx-pH', res_DEM$contrast), ]
  res_DEM.sub <- res_DEM.sub[order(res_DEM.sub$p), ]
  res_DEM.sub$text <- ifelse(res_DEM.sub$q < 0.001, '***', 
                             ifelse(res_DEM.sub$q < 0.01, '**',
                                    ifelse(res_DEM.sub$q < 0.05, '*', 'NS')))
  
  
  dat.sig <- dat.met.norm[which(rownames(dat.met.norm) %in% unique(res_DEM.sub$Index)[1:10]), ]
  dat.sig <- as.data.frame(t(dat.sig))
  dat.sig$sample <- rownames(dat.sig)
  dat.sig_df <- melt(dat.sig)
  
  dat.sig_df$Group <- gsub('^([^-]+)\\-.*$', '\\1', dat.sig_df$sample, perl = T)
  dat.sig_df$Time <- gsub('^.*\\-(\\d+)\\-.*$', '\\1', dat.sig_df$sample, perl = T)
  dat.sig_df$Time[grep('sham', dat.sig_df$Time)] <- '0'
  dat.sig_df$name <- dat.compound_DE$Compounds[match(dat.sig_df$variable, dat.compound_DE$Index)]
  num <- 1
  for(me in unique(res_DEM.sub$Index)[1:10]){
    dat.sig_df.sub <- dat.sig_df[dat.sig_df$variable == me, ]
    
    len <- max(dat.sig_df.sub$value) - min(dat.sig_df.sub$value)
    
    add_loc <- res_DEM.sub[res_DEM.sub$Index == me & res_DEM.sub$q < 0.05,]
    
    
    boxplot <- ggplot(data = dat.sig_df.sub, aes(x = Time, y = value)) +
      geom_boxplot(aes(fill = Group), outlier.size = 1) +
      labs(y = "Normalized scaled intensity", title = unique(dat.sig_df.sub$name)) +
      scale_y_continuous(limits = c(min(dat.sig_df.sub$value), min(dat.sig_df.sub$value) + len * 1.06)) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 11, hjust = .5, vjust = .5),
        strip.text.x = element_text(size = 10),
        legend.title = element_blank()
      )
    
    for(i in 1:nrow(add_loc)){
      if(grepl('18', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 2, y = min(dat.sig_df.sub$value) + len * 1.05, label = add_loc$text[i], size = 8)
      }else if(grepl('36', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 3, y = min(dat.sig_df.sub$value) + len * 1.05, label = add_loc$text[i], size = 8)
      }else  if(grepl('72', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 4, y = min(dat.sig_df.sub$value) + len * 1.05, label = add_loc$text[i], size = 8)
      }
      
    }

    ggsave(plot = boxplot, width = 4.5, height = 4, 
           filename = paste0('runtime/4.DEM.new/Treatment/Cluster',clu,'.', me, '_', num, '.boxplot.pdf'))
    num <- num + 1
    
  }
  

 
  
  
}






