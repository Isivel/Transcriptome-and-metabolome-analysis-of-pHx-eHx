
outDir <- './runtime/3.1.PLSDA'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE) 

mSet$dataSet$cls.num <- length(unique(mSet$dataSet$cls))

mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet$analSet$plsr.reg <- TRUE

mSet<-PLSDA.CV(mSet, "L",5, "Q2")

mSet<-PlotPLS.Classification(mSet, "runtime/3.1.PLSDA/pls_cv_0_", "pdf")
mSet<-PlotPLS.Imp(mSet, "runtime/3.1.PLSDA/pls_imp_0_", "pdf", type = "vip", feat.nm = "Comp. 1", feat.num = 15)
#mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)





