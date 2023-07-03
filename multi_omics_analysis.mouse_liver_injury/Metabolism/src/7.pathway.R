
#remotes::install_github("YuLab-SMU/createKEGGdb")
#library(createKEGGdb)
#create_kegg_db('mmu')

# get compound to pathway with KEGG API
#curl -g -s -S http://rest.kegg.jp/link/pathway/compound 
# get compound included in mouse
#curl -g -s -S http://rest.kegg.jp/list/compound/mmu 



source('dot_plot.R')

pathway_info <- read.delim('pathway_mouse.txt', stringsAsFactors = F, check.names = F, header = F)
colnames(pathway_info) <- c('ID', 'name')
pathway_info$ID_alias <- gsub('path:mmu', 'ko', pathway_info$ID)
pathway_info$name_alias <- gsub(' - Mus musculus \\(mouse\\)', '', pathway_info$name)


# compound2pathway <- read.delim('compound2pathway.txt', stringsAsFactors = F, check.names = F, header = F)
# colnames(compound2pathway) <- c('cpdID', 'pathway')
# compound2pathway$ko <- gsub('path:map', 'ko', compound2pathway$pathway)
# compound2pathway$cpdID <- gsub('cpd:', '', compound2pathway$cpdID)
# 
# num_compound <- length(unique(compound2pathway$cpdID))

dat_norm <- readRDS('runtime/1.data/all_sample.norm.RDS')
dat_all <- readRDS('runtime/1.data/all_sample.withAnno.RDS')

dat_all <- dat_all[dat_all$Index %in% rownames(dat_norm) & grepl('\\w',dat_all$kegg_map),]
paths <- as.character(do.call('c',sapply(dat_all$kegg_map, function(x){unlist(strsplit(x, ','))})))
paths <- paths[!grepl('-', paths)]
path_tab <- data.frame(table(paths), stringsAsFactors = F) 

##### time
outDir <- './runtime/8.pathway/Time'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

for(clu in c('1','2', '3', '4', '5', '6')){
  compound_cluster6.time <- read.delim(paste0('runtime/7.DEM/Time/Time_cluster', clu, '.tsv'),
                                       stringsAsFactors = F, check.names = F)
  compound_cluster6.time <- compound_cluster6.time[grep('\\w+',compound_cluster6.time$kegg_map, perl = T),]
  path_cluster <- as.character(do.call('c',as.list(sapply(compound_cluster6.time$kegg_map, function(x){unlist(strsplit(x, ','))}))))
  path_cluster <- path_cluster[!grepl('-', path_cluster)]
  df_path <- data.frame(table(path_cluster), stringsAsFactors = F) 
  df_path$p <- NA #p value
  df_path$RF <- NA #rich factor
  
  for(i in 1:nrow(df_path)){
    num_target_path <- path_tab$Freq[path_tab$paths %in% df_path$path_cluster[i]]
    p <- phyper(df_path$Freq[i] - 1, num_target_path, nrow(dat_all) - num_target_path, nrow(compound_cluster6.time), lower.tail = F)
    df_path$p[i] <- signif(p, digits = 4)
    df_path$RF[i] <- signif(df_path$Freq[i] / num_target_path, digits = 4)
  }
  
  df_path$name <- pathway_info$name_alias[match(df_path$path_cluster, pathway_info$ID_alias)]
  df_path <- df_path[!is.na(df_path$name),]
  
  df_path <- df_path[order(df_path$p, decreasing = F),]
  
  df_path.sub <- df_path
  if(nrow(df_path.sub) > 15) df_path.sub <- df_path[1:15,]
  
  df_path.sub$name <- factor(df_path.sub$name, levels = rev(df_path.sub$name))
  
  dot_plot <- ggplot(df_path.sub, aes(x=RF, y= name, color = p)) +
    geom_point(aes(size = Freq))
    #scale_color_gradient(low="red", high="blue") +
  
  if(max(df_path.sub$p) < 0.05){
    dot_plot <- dot_plot + 
      scale_color_gradientn(
        colors=c("red", 'grey'),
        values=rescale(c(min(df_path.sub$p), 0.05)),
        limits=c(min(df_path.sub$p), 0.05)
      )
  }else if(min(df_path.sub$p) > 0.05){
    dot_plot <- dot_plot + 
      scale_color_gradientn(
        colors=c('grey', 'blue'),
        values=rescale(c(0.05, max(df_path.sub$p))),
        limits=c(0.05, max(df_path.sub$p))
      )
  }else{
    dot_plot <- dot_plot + 
    scale_color_gradientn(
      colors=c("red", 'grey',"blue"),
      values=rescale(c(min(df_path.sub$p), 0.05,max(df_path.sub$p))),
      limits=c(min(df_path.sub$p), max(df_path.sub$p))
    )
  }
  
  dot_plot <- dot_plot + 
    labs(x = "Rich Factor", title = paste0('Cluster', clu), color = 'Pvalue', size = 'Count') +
    scale_y_discrete(labels=function(x) str_wrap(x,width=50)) +
    guides(
      size = guide_legend(order = 2),
      color = guide_colorbar(order = 1)
    ) +
    theme_bw() +
    theme(axis.text.y = element_text( size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = .5, size = 15))
  
  ggsave(dot_plot, width = 7.5, height = 6.5,
         filename = paste0('runtime/8.pathway/Time/Time.Cluster', clu, '.keggEnrich.pdf'))
  
  write.table(df_path, col.names = T, row.names = F, sep = "\t", quote = F,
              file = paste0('runtime/8.pathway/Time/Time.Cluster', clu, '.keggEnrich.tsv'))
  
}




##### treatment
outDir <- './runtime/8.pathway/Treatment'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

for(clu in c('1','2', '3')){
  compound_cluster6.time <- read.delim(paste0('runtime/7.DEM/Treatment/treatment_Cluster', clu, '.tsv'),
                                       stringsAsFactors = F, check.names = F)
  compound_cluster6.time <- compound_cluster6.time[grep('\\w+',compound_cluster6.time$kegg_map, perl = T),]
  path_cluster <- as.character(do.call('c',as.list(sapply(compound_cluster6.time$kegg_map, function(x){unlist(strsplit(x, ','))}))))
  path_cluster <- path_cluster[!grepl('-', path_cluster)]
  df_path <- data.frame(table(path_cluster), stringsAsFactors = F) 
  df_path$p <- NA #p value
  df_path$RF <- NA #rich factor
  
  for(i in 1:nrow(df_path)){
    num_target_path <- path_tab$Freq[path_tab$paths %in% df_path$path_cluster[i]]
    p <- phyper(df_path$Freq[i] - 1, num_target_path, nrow(dat_all) - num_target_path, nrow(compound_cluster6.time), lower.tail = F)
    df_path$p[i] <- signif(p, digits = 4)
    df_path$RF[i] <- signif(df_path$Freq[i] / num_target_path, digits = 4)
  }
  
  df_path$name <- pathway_info$name_alias[match(df_path$path_cluster, pathway_info$ID_alias)]
  df_path <- df_path[!is.na(df_path$name),]
  
  df_path <- df_path[order(df_path$p, decreasing = F),]
  
  df_path.sub <- df_path
  if(nrow(df_path.sub) > 15) df_path.sub <- df_path[1:15,]
  
  df_path.sub$name <- factor(df_path.sub$name, levels = rev(df_path.sub$name))
  
  dot_plot <- ggplot(df_path.sub, aes(x=RF, y= name, color = p)) +
    geom_point(aes(size = Freq))
  #scale_color_gradient(low="red", high="blue") +
  
  if(max(df_path.sub$p) < 0.05){
    dot_plot <- dot_plot + 
      scale_color_gradientn(
        colors=c("red", 'grey'),
        values=rescale(c(min(df_path.sub$p), 0.05)),
        limits=c(min(df_path.sub$p), 0.05)
      )
  }else if(min(df_path.sub$p) > 0.05){
    dot_plot <- dot_plot + 
      scale_color_gradientn(
        colors=c('grey', 'blue'),
        values=rescale(c(0.05, max(df_path.sub$p))),
        limits=c(0.05, max(df_path.sub$p))
      )
  }else{
    dot_plot <- dot_plot + 
      scale_color_gradientn(
        colors=c("red", 'grey',"blue"),
        values=rescale(c(min(df_path.sub$p), 0.05,max(df_path.sub$p))),
        limits=c(min(df_path.sub$p), max(df_path.sub$p))
      )
  }
  
  dot_plot <- dot_plot + 
    labs(x = "Rich Factor", title = paste0('Cluster', clu), color = 'Pvalue', size = 'Count') +
    scale_y_discrete(labels=function(x) str_wrap(x,width=50)) +
    guides(
      size = guide_legend(order = 2),
      color = guide_colorbar(order = 1)
    ) +
    theme_bw() +
    theme(axis.text.y = element_text( size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = .5, size = 15))
  
  ggsave(dot_plot, width = 7.5, height = 6.5,
         filename = paste0('runtime/8.pathway/Treatment/Treatment.Cluster', clu, '.keggEnrich.pdf'))
  
  write.table(df_path, col.names = T, row.names = F, sep = "\t", quote = F,
              file = paste0('runtime/8.pathway/Treatment/Treatment.Cluster', clu, '.keggEnrich.tsv'))
  
}


##### enrichment for combined all treatment
files <- list.files(path = 'runtime/7.DEM/Treatment', pattern = '^Treatment_cluster\\d.tsv', full.names = T)
compound_cluster6.time <- NULL
for(f in files){
  dat <- read.delim(f, stringsAsFactors = F, check.names = F)
  compound_cluster6.time <- rbind(compound_cluster6.time, dat)
}
compound_cluster6.time <- compound_cluster6.time[grep('\\w+',compound_cluster6.time$kegg_map, perl = T),]
path_cluster <- as.character(do.call('c',as.list(sapply(compound_cluster6.time$kegg_map, function(x){unlist(strsplit(x, ','))}))))
path_cluster <- path_cluster[!grepl('-', path_cluster)]
df_path <- data.frame(table(path_cluster), stringsAsFactors = F) 
df_path$p <- NA #p value
df_path$RF <- NA #rich factor

for(i in 1:nrow(df_path)){
  num_target_path <- path_tab$Freq[path_tab$paths %in% df_path$path_cluster[i]]
  p <- phyper(df_path$Freq[i] - 1, num_target_path, nrow(dat_all) - num_target_path, nrow(compound_cluster6.time), lower.tail = F)
  df_path$p[i] <- signif(p, digits = 4)
  df_path$RF[i] <- signif(df_path$Freq[i] / num_target_path, digits = 4)
}

df_path$name <- pathway_info$name_alias[match(df_path$path_cluster, pathway_info$ID_alias)]
df_path <- df_path[!is.na(df_path$name),]

df_path <- df_path[order(df_path$p, decreasing = F),]

df_path.sub <- df_path
if(nrow(df_path.sub) > 15) df_path.sub <- df_path[1:15,]

df_path.sub$name <- factor(df_path.sub$name, levels = rev(df_path.sub$name))

dot_plot <- ggplot(df_path.sub, aes(x=RF, y= name, color = p)) +
  geom_point(aes(size = Freq))
#scale_color_gradient(low="red", high="blue") +

if(max(df_path.sub$p) < 0.05){
  dot_plot <- dot_plot + 
    scale_color_gradientn(
      colors=c("red", 'grey'),
      values=rescale(c(min(df_path.sub$p), 0.05)),
      limits=c(min(df_path.sub$p), 0.05)
    )
}else if(min(df_path.sub$p) > 0.05){
  dot_plot <- dot_plot + 
    scale_color_gradientn(
      colors=c('grey', 'blue'),
      values=rescale(c(0.05, max(df_path.sub$p))),
      limits=c(0.05, max(df_path.sub$p))
    )
}else{
  dot_plot <- dot_plot + 
    scale_color_gradientn(
      colors=c("red", 'grey',"blue"),
      values=scales::rescale(c(min(df_path.sub$p), 0.05,max(df_path.sub$p))),
      limits=c(min(df_path.sub$p), max(df_path.sub$p))
    )
}

dot_plot <- dot_plot + 
  labs(x = "Rich Factor", title = "Differential metabolites in treatment", color = 'Pvalue', size = 'Count') +
  scale_y_discrete(labels=function(x) str_wrap(x,width=50)) +
  guides(
    size = guide_legend(order = 2),
    color = guide_colorbar(order = 1)
  ) +
  theme_bw() +
  theme(axis.text.y = element_text( size = 11),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = .5, size = 15))

ggsave(dot_plot, width = 7.5, height = 6.5,
        filename = paste0('runtime/7.DEM/Treatment/Treatment.keggEnrich.pdf'))

write.table(df_path, col.names = T, row.names = F, sep = "\t", quote = F,
            file = paste0('runtime/7.DEM/Treatment/Treatment.keggEnrich.tsv'))


write.table(compound_cluster6.time, col.names = T, row.names = T, sep = "\t", quote = F,
            file = 'runtime/7.DEM/Treatment/Treatment_DEM.tsv')



