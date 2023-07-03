
library(ggplot2)

outDir <- './runtime/4.DEM.new'
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

res_DEM <- NULL
for(me in rownames(dat.met.norm.twA)){
  dat_df <- df_calDiff
  dat_df$expr <- as.numeric(dat.met.norm.twA[me,])
  dat_df$Group <- factor(dat_df$Group, levels = c('pHx', 'eHx'))
  times1 <- c('0', '18', '36', '72')
  for(ti in times1){
    if(ti != '0'){
      res_tmp <- t.test(as.numeric(dat_df$expr[dat_df$Group == 'eHx' & dat_df$Time == ti]),
                        as.numeric(dat_df$expr[dat_df$Group == 'pHx' & dat_df$Time == ti]))
      res_DEM <- rbind(res_DEM, c(me, paste0('eHx-pHx:', ti),
                                  signif(c(res_tmp$statistic, 
                                           as.numeric(res_tmp$estimate), 
                                           res_tmp$p.value), digits = 4)))
    }

    index <- which(times1 %in% ti)
    times2 <- times1[-c(1:index)]
    if(length(times2) == 0) next()
    for(tii in times2){
      if(ti == tii) next()
      res_pH <- t.test(as.numeric(dat_df$expr[dat_df$Group == 'pHx' & dat_df$Time == ti]),
                        as.numeric(dat_df$expr[dat_df$Group == 'pHx' & dat_df$Time == tii]))
      res_DEM <- rbind(res_DEM, c(me, paste0(ti, '-', tii, ':pHx'),
                                  signif(c(res_pH$statistic, 
                                           as.numeric(res_pH$estimate), 
                                           res_pH$p.value), digits = 4)))
      
      
      res_eHx <- t.test(as.numeric(dat_df$expr[dat_df$Group == 'eHx' & dat_df$Time == ti]),
                       as.numeric(dat_df$expr[dat_df$Group == 'eHx' & dat_df$Time == tii]))
      res_DEM <- rbind(res_DEM, c(me, paste0(ti, '-', tii, ':eHx'),
                                  signif(c(res_eHx$statistic, 
                                           as.numeric(res_eHx$estimate), 
                                           res_eHx$p.value), digits = 4)))
    }
    
  }
  
  
  
  #model<-aov(expr ~ Group * Time, data=dat_df)
  #res.anova <- summary(model)
  #res.out <- LSD.test(model,"Time",p.adj="hommel",group = T, console = F)
  #res.out <- LSD.test(model,"Group",p.adj="hommel",group = T, console = F)
  #if(length(unique(res.out$groups$groups)) > 1 & res.anova[[1]][1,5] < 0.01){ 
  #  sig_me <- rbind(sig_me, c(as.numeric(res.out$groups[c('0','18','36','72'), 'expr']), res.anova[[1]][1,5]))
  #  compound <- c(compound, me)
  #}
  
  #res_group <- TukeyHSD(model, which = 'Group')
  #res_time <- TukeyHSD(model, which = 'Time')
  

  
  #res_sig <- rbind(res_sig, c(res.anova[[1]][1,5], res.anova[[1]][2,5]))
}

res_DEM <- as.data.frame(res_DEM, stringsAsFactors = F)
colnames(res_DEM) <- c('Index','contrast', 't', 'mean.g1', 'mean.g2', 'p')
storage.mode(res_DEM$t) <- 'numeric'
storage.mode(res_DEM$mean.g1) <- 'numeric'
storage.mode(res_DEM$mean.g2) <- 'numeric'
storage.mode(res_DEM$p) <- 'numeric'
res_DEM <- res_DEM[order(res_DEM$contrast), ]

res_DEM$q <- c(sapply(unique(res_DEM$contrast), function(x){
  as.numeric(p.adjust(res_DEM$p[which(res_DEM$contrast == x)], method = 'fdr'))
}))

dat_anno <- readRDS('runtime/1.data/all_sample.withAnno.RDS')
res_DEM$Compound <- dat_anno$Compounds[match(res_DEM$Index, dat_anno$Index)]
res_DEM$Class.I <- dat_anno$Class.I[match(res_DEM$Index, dat_anno$Index)]
res_DEM$Class.II <- dat_anno$Class.II[match(res_DEM$Index, dat_anno$Index)]

save(res_DEM, file = 'runtime/4.DEM.new/res_DEM.rda')

write.table(res_DEM, col.names = T, row.names = F, sep = "\t", quote = F,
            file = "runtime/4.DEM.new/res_DEM.all.tsv")

