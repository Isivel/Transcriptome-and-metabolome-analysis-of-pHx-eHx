
library(ropls)
library(ggplot2)
library(MetaboAnalystR)
library(agricolae)
library(pls)

outDir <- './runtime/2.PCA'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE) 

##### processed MetaboAnalystR
mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "runtime/1.data/All_rawdat.csv", "colu", "disc")

dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')

mSet$dataSet$norm <- t(dat.met.norm)


#### PCA
mSet <- PCA.Anal(mSet)
mSet <- PlotPCAPairSummary(mSet, "runtime/2.PCA/pca_pair_0_", "pdf", 72, width=NA, 5)
mSet <- PlotPCAScree(mSet, "runtime/2.PCA/pca_scree_0_", "pdf", 72, width=NA, 5)
mSet <- PlotPCA2DScore(mSet, "runtime/2.PCA/pca_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,1,0)
mSet <- PlotPCALoading(mSet, "runtime/2.PCA/pca_loading_0_", "pdf", 72, width=NA, 1,2);
mSet <- PlotPCABiplot(mSet, "runtime/2.PCA/pca_biplot_0_", "pdf", 72, width=NA, 1,2)
mSet <- PlotPCA3DScoreImg(mSet, "runtime/2.PCA/pca_score3d_0_", "pdf", 72, width=NA, 1,2,3, 40)

save(mSet, file = 'runtime/2.PCA/mSet.Rdata')

