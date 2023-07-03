
library(ropls)
library(ggplot2)
library(MetaboAnalystR)
library(agricolae)
library(pls)

outDir <- './runtime/6.VIP_score'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE) 

##### Time series in eHx
dat <- read.csv('runtime/1.data/All_rawdat.csv', stringsAsFactors = F, check.names = F)
dat.sub <- dat[, !grepl('pHx', colnames(dat))]
write.table(dat.sub, col.names = T, row.names = F, sep = ",", quote = F,
            file = "runtime/6.VIP_score/Time.eHx.rawData.csv")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "runtime/6.VIP_score/Time.eHx.rawData.csv", "colu", "disc")

dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')
dat.met.norm.sub <- dat.met.norm[, !grepl('pHx', colnames(dat.met.norm))]

mSet$dataSet$norm <- t(dat.met.norm.sub)

mSet <- PLSR.Anal(mSet, reg=TRUE)
mSet$analSet$plsr.reg <- TRUE
mSet$dataSet$cls.num <- length(unique(mSet$dataSet$cls))
mSet<-PLSDA.CV(mSet, "L",5, "Q2")

outDir <- './runtime/6.VIP_score/Time.eHx'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

mSet<-PlotPLS.Classification(mSet, "runtime/6.VIP_score/Time.eHx/pls_cv_0_", "pdf")
mSet<-PlotPLS.Imp(mSet, "runtime/6.VIP_score/Time.eHx/pls_imp_0_", "pdf", type = "vip", feat.nm = "Comp. 1", feat.num = 15)

saveRDS(mSet$analSet$plsda$vip.mat, file = 'runtime/6.VIP_score/Time.eHx/vip.mat.rds')
save(mSet, file = 'runtime/6.VIP_score/Time.eHx/mSet.rda')

rm(mSet)

#### Time series for pHx
dat.sub <- dat[, !grepl('eHx', colnames(dat))]
write.table(dat.sub, col.names = T, row.names = F, sep = ",", quote = F,
            file = "runtime/6.VIP_score/Time.pHx.rawData.csv")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "runtime/6.VIP_score/Time.pHx.rawData.csv", "colu", "disc")

dat.met.norm.sub <- dat.met.norm[, !grepl('eHx', colnames(dat.met.norm))]

mSet$dataSet$norm <- t(dat.met.norm.sub)

mSet <- PLSR.Anal(mSet, reg=TRUE)
mSet$analSet$plsr.reg <- TRUE
mSet$dataSet$cls.num <- length(unique(mSet$dataSet$cls))
mSet<-PLSDA.CV(mSet, "L",5, "Q2")

outDir <- './runtime/6.VIP_score/Time.pHx'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

mSet<-PlotPLS.Classification(mSet, "runtime/6.VIP_score/Time.pHx/pls_cv_0_", "pdf")
mSet<-PlotPLS.Imp(mSet, "runtime/6.VIP_score/Time.pHx/pls_imp_0_", "pdf", type = "vip", feat.nm = "Comp. 1", feat.num = 15)

saveRDS(mSet$analSet$plsda$vip.mat, file = 'runtime/6.VIP_score/Time.pHx/vip.mat.rds')
save(mSet, file = 'runtime/6.VIP_score/Time.pHx/mSet.rda')


##### Condition in time 18h
dat.sub <- dat[, c(1, grep('18', colnames(dat)))]
dat.sub[1,grep('pHx', colnames(dat.sub))] <- 'pHx'
dat.sub[1,grep('eHx', colnames(dat.sub))] <- 'eHx'
write.table(dat.sub, col.names = T, row.names = F, sep = ",", quote = F,
            file = "runtime/6.VIP_score/Treatment.18.rawData.csv")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "runtime/6.VIP_score/Treatment.18.rawData.csv", "colu", "disc")

dat.met.norm.sub <- dat.met.norm[, grep('18', colnames(dat.met.norm))]

mSet$dataSet$norm <- t(dat.met.norm.sub)

mSet <- PLSR.Anal(mSet, reg=TRUE)
mSet$analSet$plsr.reg <- TRUE
mSet$dataSet$cls.num <- length(unique(mSet$dataSet$cls))
mSet<-PLSDA.CV(mSet, "L",5, "Q2")

outDir <- './runtime/6.VIP_score/Treatment.18'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

mSet<-PlotPLS.Classification(mSet, "runtime/6.VIP_score/Treatment.18/pls_cv_0_", "pdf")
mSet<-PlotPLS.Imp(mSet, "runtime/6.VIP_score/Treatment.18/pls_imp_0_", "pdf", type = "vip", feat.nm = "Comp. 1", feat.num = 15)

saveRDS(mSet$analSet$plsda$vip.mat, file = 'runtime/6.VIP_score/Treatment.18/vip.mat.rds')
save(mSet, file = 'runtime/6.VIP_score/Treatment.18/mSet.rda')


##### Condition in time 36h
rm(mSet)
dat.sub <- dat[, c(1, grep('36', colnames(dat)))]
dat.sub[1,grep('pHx', colnames(dat.sub))] <- 'pHx'
dat.sub[1,grep('eHx', colnames(dat.sub))] <- 'eHx'
write.table(dat.sub, col.names = T, row.names = F, sep = ",", quote = F,
            file = "runtime/6.VIP_score/Treatment.36.rawData.csv")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "runtime/6.VIP_score/Treatment.36.rawData.csv", "colu", "disc")

dat.met.norm.sub <- dat.met.norm[, grep('36', colnames(dat.met.norm))]

mSet$dataSet$norm <- t(dat.met.norm.sub)

mSet <- PLSR.Anal(mSet, reg=TRUE)
mSet$analSet$plsr.reg <- TRUE
mSet$dataSet$cls.num <- length(unique(mSet$dataSet$cls))
mSet<-PLSDA.CV(mSet, "L",5, "Q2")

outDir <- './runtime/6.VIP_score/Treatment.36'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

mSet<-PlotPLS.Classification(mSet, "runtime/6.VIP_score/Treatment.36/pls_cv_0_", "pdf")
mSet<-PlotPLS.Imp(mSet, "runtime/6.VIP_score/Treatment.36/pls_imp_0_", "pdf", type = "vip", feat.nm = "Comp. 1", feat.num = 15)

saveRDS(mSet$analSet$plsda$vip.mat, file = 'runtime/6.VIP_score/Treatment.36/vip.mat.rds')
save(mSet, file = 'runtime/6.VIP_score/Treatment.36/mSet.rda')


##### Condition in time 72h
rm(mSet)
dat.sub <- dat[, c(1, grep('72', colnames(dat)))]
dat.sub[1,grep('pHx', colnames(dat.sub))] <- 'pHx'
dat.sub[1,grep('eHx', colnames(dat.sub))] <- 'eHx'
write.table(dat.sub, col.names = T, row.names = F, sep = ",", quote = F,
            file = "runtime/6.VIP_score/Treatment.72.rawData.csv")

mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "runtime/6.VIP_score/Treatment.72.rawData.csv", "colu", "disc")

dat.met.norm.sub <- dat.met.norm[, grep('72', colnames(dat.met.norm))]

mSet$dataSet$norm <- t(dat.met.norm.sub)

mSet <- PLSR.Anal(mSet, reg=TRUE)
mSet$analSet$plsr.reg <- TRUE
mSet$dataSet$cls.num <- length(unique(mSet$dataSet$cls))
mSet<-PLSDA.CV(mSet, "L",5, "Q2")

outDir <- './runtime/6.VIP_score/Treatment.72'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

mSet<-PlotPLS.Classification(mSet, "runtime/6.VIP_score/Treatment.72/pls_cv_0_", "pdf")
mSet<-PlotPLS.Imp(mSet, "runtime/6.VIP_score/Treatment.72/pls_imp_0_", "pdf", type = "vip", feat.nm = "Comp. 1", feat.num = 15)

saveRDS(mSet$analSet$plsda$vip.mat, file = 'runtime/6.VIP_score/Treatment.72/vip.mat.rds')
save(mSet, file = 'runtime/6.VIP_score/Treatment.72/mSet.rda')










