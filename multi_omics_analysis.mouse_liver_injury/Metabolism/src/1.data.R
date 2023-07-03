

library(ropls)
library(ggplot2)

### as pipeline https://www.intechopen.com/chapters/52527

outDir <- './runtime/1.data'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE) 

dat_anno <- openxlsx::read.xlsx('ALL_sample_data.xlsx', check.names = F)
dat.met <- dat_anno[, c(1,12:51)]
colnames(dat_anno) <- gsub('^pH-', 'pHx-', colnames(dat_anno), perl = T)
colnames(dat.met) <- gsub('^pH-', 'pHx-', colnames(dat.met), perl = T)

rownames(dat.met) <- dat.met$Index
dat.met$Index <- NULL


#### remove features with low repeatability in QC samples
cal_cv <- function(x){
  x <- as.numeric(x)
  return(sd(x) / mean(x))
}

res_cv <- apply(dat.met[,36:40], 1, cal_cv)

dat.met <- dat.met[res_cv < 0.3, 1:35]



### prepare for MetaboAnalystR
group <- gsub('^([^-]+)-.*$', '\\1', colnames(dat.met), perl = T)
time <- gsub('^.*-(\\d+)-.*$', '\\1', colnames(dat.met), perl = T)
time[grep('sham', time)] <- '0'
dat_tmp <- rbind(time, dat.met)
rownames(dat_tmp)[1] <- "Label"
write.csv(data.frame(Sample = rownames(dat_tmp), dat_tmp), row.names = F,
          file = 'runtime/1.data/All_rawdat.csv')


paretoscale <- function(z) {
  z <- log2(z + 1)
  rowmean <- apply(z, 1, mean) # row means
  rowsd <- apply(z, 1, sd) # row standard deviation
  rowsqrtsd <- sqrt(rowsd) # sqrt of sd
  rv <- sweep(z, 1, rowmean,"-") # mean center
  rv <- sweep(rv, 1, rowsqrtsd, "/") # divide by sqrtsd
  return(rv)}

dat.met.norm <- paretoscale(dat.met)


saveRDS(dat_anno, file = 'runtime/1.data/all_sample.withAnno.RDS')
saveRDS(dat.met, file = 'runtime/1.data/all_sample.onlyExpr.RDS')
saveRDS(dat.met.norm, file = 'runtime/1.data/all_sample.norm.RDS')






