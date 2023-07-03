
library(ropls)
library(ggplot2)
library(MetaboAnalystR)
library(agricolae)
library(pls)

outDir <- './runtime/3.OPLSDA'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE) 

load('runtime/2.PCA/mSet.Rdata')

#### OPLSDA
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet<-PlotOPLS2DScore(mSet, "runtime/3.OPLSDA/opls_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotOPLS.Splot(mSet, "runtime/3.OPLSDA/opls_splot_0_", "all", "pdf", 72, width=NA);
mSet<-PlotOPLS.MDL(mSet, "runtime/3.OPLSDA/opls_mdl_0_", "pdf", 72, width=NA)
#mSet<-PlotOPLS.Imp(mSet, "runtime/3.OPLSDA/opls_imp_0_", "pdf", 72, width=NA, "vip", "tscore", 15,FALSE)
mSet$analSet$opls.reg <- TRUE
mSet<-OPLSDA.Permut(mSet, 1000)
mSet<-PlotOPLS.Permutation(mSet, "runtime/3.OPLSDA/opls_perm_1_", "pdf", 72, width=NA)


save(mSet, file = 'runtime/3.OPLSDA/mSet.Rdata')

load('runtime/3.OPLSDA/mSet.Rdata')


