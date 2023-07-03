

library(ggplot2)
library(reshape2)

outDir <- './runtime/7.DEM/Time'
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
df_calDiff$Group[1:5] <- 'pHx'
df_calDiff$Group[6:10] <- 'eHx'


load('runtime/4.DEM.new/res_DEM.rda')

vip1 <- readRDS('runtime/6.VIP_score/Time.eHx/vip.mat.rds')
vip2 <- readRDS('runtime/6.VIP_score/Time.pHx/vip.mat.rds')

compounds <- intersect(unique(res_DEM$Index[which(grepl('^\\d', res_DEM$contrast, perl = T) & res_DEM$q < 0.05)]),
                       c(rownames(vip1)[which(rowMeans(vip1[, 1:2]) > 1.5)], rownames(vip2)[which(rowMeans(vip2[, 1:2]) > 1.5)]))

plotDir <- 'runtime/7.DEM/Time'
## plot heatmap
annotation_col = df_calDiff[!grepl('sham.*_2',df_calDiff$sample, perl = T), ]
annotation_col$Group[grepl('sham', annotation_col$sample)] <- 'sham'
rownames(annotation_col) <- annotation_col$sample
annotation_col$sample <- NULL
ann_colors = list(
  Group = c("sham"="steelblue","eHx"="red", "pHx"="#FF7F00"),
  Time = c('0' = 'navy', '18' = '#ffcccb', '36' = 'salmon', '72' = 'darkred')
)
sample_order<-order(annotation_col$Group, annotation_col$Time)
pdf("runtime/7.DEM/Time/heatmap.Time.pdf", width = 7.5, height = 7.5)
out_Time <- pheatmap::pheatmap(dat.met.norm[compounds,],
                               clustering_method = "ward.D2",
                               clustering_distance_rows = 'correlation',
                               show_rownames =F,show_colnames =F,fontsize_row=7,
                               annotation_col = annotation_col,cluster_cols = F,scale = 'row',
                               annotation_colors = ann_colors,color = colorRampPalette(c("blue", "white", "firebrick3"))(100)) 
dev.off()

save(out_Time, file = "runtime/7.DEM/Time/pheatmap_Time.Rdata")


#############################################################
## extract metabolist with cluster
all_data <- readRDS('runtime/1.data/all_sample.withAnno.RDS')

meta_tab.Time <- sort(cutree(out_Time$tree_row, k = 6))


cluster_tab.time <- sort(cutree(out_Time$tree_row, k = 6))
annotation_row <- data.frame(ID = compounds, stringsAsFactors = F)
annotation_row$Cluster <- cluster_tab.time[match(annotation_row$ID, names(cluster_tab.time))]
annotation_row$Cluster <- paste0('Cluster', annotation_row$Cluster)
rownames(annotation_row) <- annotation_row$ID
annotation_row <- annotation_row[annotation_row$Cluster != 'ClusterNA',]
tmp_tab <- annotation_row
tmp_tab$Compound <- res_DEM$Compound[match(tmp_tab$ID, res_DEM$Index)]
tmp_tab$Class.I <- res_DEM$Class.I[match(tmp_tab$ID, res_DEM$Index)]
tmp_tab$Class.II <- res_DEM$Class.II[match(tmp_tab$ID, res_DEM$Index)]
write.table(tmp_tab, col.names = T, row.names = F, sep = "\t", quote = F,
            file = "runtime/7.DEM/Time/heatmap.time.cluster.txt")
annotation_row$ID <- NULL
require(randomcoloR)
n <- 6
palette <- distinctColorPalette(n)
ann_colors = list(
  Group = c("sham"="steelblue","eHx"="red", "pHx"="#FF7F00"),
  Time = c('0' = 'navy', '18' = '#ffcccb', '36' = 'salmon', '72' = 'darkred'),
  Cluster = c('Cluster1' = palette[1], 'Cluster2' = palette[2], 'Cluster3' = palette[3], 'Cluster4' = palette[4],
              'Cluster5' = palette[5], 'Cluster6' = palette[6])
)

sample_order<-order(annotation_col$Group, annotation_col$Time)
pdf("runtime/7.DEM/Time/heatmap.Time.withCluster.pdf", width = 7.5, height = 7.5)
out_Time <- pheatmap::pheatmap(dat.met.norm[compounds, ],
                               clustering_method = "ward.D2",
                               clustering_distance_rows = 'correlation',
                               show_rownames =F,show_colnames =F,fontsize_row=7,
                               annotation_col = annotation_col,annotation_row = annotation_row,
                               cluster_cols = F,scale = 'row',
                               annotation_colors = ann_colors,color = colorRampPalette(c("blue", "white", "firebrick3"))(100)) 
dev.off()




line_plot.list <- NULL
pie_plot.list <- NULL
for(clu in 1:6){
  compound.list <- all_data$Compounds[which(all_data$Index %in% names(meta_tab.Time)[meta_tab.Time == clu])]
  write.table(compound.list, row.names = F, col.names = F, quote = F, file = paste0('runtime/7.DEM/Time/compound_Time_cluster', clu,'.tsv'))
  
  compound_cluster <- all_data[which(all_data$Index %in% names(meta_tab.Time)[meta_tab.Time == clu]),]
  write.table(compound_cluster, row.names = T, col.names = T, quote = F, sep = "\t",
              file = paste0('runtime/7.DEM/Time/Time_cluster', clu,'.tsv'))
  df_plot.pie <- as.data.frame(table(compound_cluster$Class.I), stringsAsFactors = F)
  df_plot.pie$lab <- paste0(df_plot.pie$Var1, ' (', df_plot.pie$Freq, ')')
  df_plot.pie$lab <- factor(df_plot.pie$lab, levels = df_plot.pie$lab[order(df_plot.pie$Freq, decreasing = T)])
  ## pie
  pie_plot.list[[paste0('Cluster', clu)]] <- ggplot(data = df_plot.pie, aes(x="", y=Freq, fill=lab)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    labs(title = paste0('Cluster', clu)) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = .5, vjust = .5),
      legend.title = element_blank()
    )
  ggsave(plot =  pie_plot.list[[paste0('Cluster', clu)]], width = 9, height = 6, filename = paste0('runtime/7.DEM/Time/compound_cluster', clu,'.pie.pdf'))
  ## trend line
  df_tmp <- as.data.frame(t(dat.met.norm.twA[names(meta_tab.Time)[which(meta_tab.Time == clu)],]))
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
  cols[grep('pHxCol', levels(df_plot.line$ID))] <- 'steelblue'
  cols[grep('eHxCol', levels(df_plot.line$ID))] <- 'salmon'
  
  df_plot.line$al <- 0.6
  df_plot.line$al[grepl('Col', df_plot.line$ID)] <- 1
  
  df_plot.line$Time <- as.numeric(df_plot.line$Time)
  
  line_plot.list[[paste0('Cluster', clu)]] <- ggplot(data = df_plot.line, aes(x = Time, y = value, color = ID)) +
    geom_line(aes(alpha = al)) +
    scale_color_manual(values = cols) +
    scale_x_continuous(breaks = c(0, 18, 36, 72), labels = c('0','18','36','72')) +
    labs(y = 'Normalized scaled intensity', title = paste0('Cluster', clu)) +
    theme(
      plot.title = element_text(hjust = .5, vjust = .5),
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
  
  ggsave(plot = line_plot.list[[paste0('Cluster', clu)]], width = 6.5, height = 6, filename = paste0('runtime/7.DEM/Time/compound_cluster', clu,'_linePlot.pdf'))
  
}

require(ggpubr)
pdf("runtime/7.DEM/Time/DEM.line.pdf",width = 15,height = 9, onefile = T)
ggarrange(plotlist = line_plot.list,ncol = 3, nrow = 2)
dev.off()

pdf("runtime/7.DEM/Time/DEM.pie.pdf",width = 15,height = 9, onefile = T)
ggarrange(plotlist = pie_plot.list,ncol = 3, nrow = 2)
dev.off()



#### boxplot
dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')
plotDir <- './runtime/7.DEM/Time/boxplot'
if(!dir.exists(plotDir)) dir.create(plotDir, recursive = TRUE)

for(clu in c('1', '2', '3', '4', '5', '6')){
  dat.compound_DE <- all_data[which(all_data$Index %in% names(meta_tab.Time)[meta_tab.Time == clu]),]
  res_DEM.sub <- res_DEM[res_DEM$Index %in% dat.compound_DE$Index & grepl('^\\d', res_DEM$contrast, perl = T), ]
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
    
    #min(dat.sig_df.sub$value) <- min(dat.sig_df.sub$value)
    len <- max(dat.sig_df.sub$value) - min(dat.sig_df.sub$value)
    
    add_loc <- res_DEM.sub[res_DEM.sub$Index == me & res_DEM.sub$q < 0.05,]
    
    
    boxplot <- ggplot(data = dat.sig_df.sub, aes(x = Time, y = value)) +
      geom_boxplot(aes(fill = Group), outlier.size = 1) +
      labs(y = "Normalized scaled intensity", title = unique(dat.sig_df.sub$name)) +
      #scale_y_continuous(limits = c(min(dat.sig_df.sub$value), min(dat.sig_df.sub$value) + len * 1.06)) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 11, hjust = .5, vjust = .5),
        strip.text.x = element_text(size = 10),
        legend.title = element_blank()
      )
    
    
    for(i in 1:nrow(add_loc)){
      Y1 <- min(dat.sig_df.sub$value) + len * (1.03 + (i - 1) * 0.03)
      Y2 <- min(dat.sig_df.sub$value) + len * (1.02 + (i - 1) * 0.03)
      if(grepl('0-18:eHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1 + (1.75 - 1)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1, xend = 1.75, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1, xend = 1, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1.75, xend = 1.75, 
                                  y = Y2, 
                                  yend = Y1))
        
      }else if(grepl('0-18:pHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1 + (2.25 - 1)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1, xend = 2.25, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1, xend = 1, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 2.25, xend = 2.25, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('0-36:eHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1 + (2.75 - 1)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1, xend = 2.75, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1, xend = 1, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 2.75, xend = 2.75, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('0-36:pHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1 + (3.25 - 1)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1, xend = 3.25, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1, xend = 1, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 3.25, xend = 3.25, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('0-72:eHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1 + (3.75 - 1)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1, xend = 3.75, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1, xend = 1, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 3.75, xend = 3.75, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('0-72:pHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1 + (4.25 - 1)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1, xend = 4.25, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1, xend = 1, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 4.25, xend = 4.25, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('18-36:eHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1.75 + (2.75 - 1.75)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1.75, xend = 2.75, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1.75, xend = 1.75, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 2.75, xend = 2.75, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('18-36:pHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 2.25 + (3.25 - 2.25)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 2.25, xend = 3.25, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 2.25, xend = 2.25, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 3.25, xend = 3.25, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('18-72:eHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 1.75 + (3.75 - 1.75)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 1.75, xend = 3.75, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 1.75, xend = 1.75, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 3.75, xend = 3.75, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('18-72:pHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 2.25 + (4.25 - 2.25)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 2.25, xend = 4.25, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 2.25, xend = 2.25, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 4.25, xend = 4.25, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('36-72:eHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 2.75 + (3.75 - 2.75)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 2.75, xend = 3.75, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 2.75, xend = 2.75, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 3.75, xend = 3.75, 
                                  y = Y2, 
                                  yend = Y1))
      }else if(grepl('36-72:pHx', add_loc$contrast[i])){
        boxplot <- boxplot + annotate('text', x = 3.25 + (4.25 - 3.25)/2, 
                                      y = Y1, 
                                      label = add_loc$text[i], size = 6) +
          geom_segment(aes_string(x = 3.25, xend = 4.25, 
                                  y = Y1, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 3.25, xend = 3.25, 
                                  y = Y2, 
                                  yend = Y1)) +
          geom_segment(aes_string(x = 4.25, xend = 4.25, 
                                  y = Y2, 
                                  yend = Y1))
      }
      
    }
    
    ggsave(plot = boxplot, width = 4.5, height = 4.5, 
           filename = paste0('runtime/7.DEM/Time/boxplot/Cluster',clu,'.', me, '_', num, '.boxplot.pdf'))
    num <- num + 1
    
    
  }
  
}






