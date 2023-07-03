
library(ggvenn)

outDir <- './runtime/8.venn'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
# 
# load('runtime/4.DEM.new/res_DEM.rda')
# 
# vip1 <- readRDS('runtime/6.VIP_score/Time.eHx/vip.mat.rds')
# vip2 <- readRDS('runtime/6.VIP_score/Time.pHx/vip.mat.rds')
# 
# 
# contrasts_target.pHx <- c('0-18:pHx', '0-36:pHx', '0-72:pHx')
# contrasts_target.eHx <- c('0-18:eHx', '0-36:eHx', '0-72:eHx')
# 
# res_DEM.pHx <- res_DEM[res_DEM$q < 0.05 & 
#                          res_DEM$Index %in% rownames(vip2)[which(rowMeans(vip2[, 1:2]) > 1.5)] &
#                          res_DEM$contrast %in% contrasts_target.pHx, c('Index','contrast')]
# res_DEM.eHx <- res_DEM[res_DEM$q < 0.05 & 
#                          res_DEM$Index %in% rownames(vip2)[which(rowMeans(vip1[, 1:2]) > 1.5)] &
#                          res_DEM$contrast %in% contrasts_target.eHx, c('Index','contrast')]
# 
# 
# df_plot <- rbind(res_DEM.eHx, res_DEM.pHx)
# 
# split_tibble <- function(tibble, column = 'col') {
#   tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
# }
# x <- split_tibble(df_plot, 'contrast')
# 
# mycolor <- c("blue", "yellow", "steelblue", "red", 'salmon', "green")
# 
# venn_plot <- ggvenn(
#   x, 
#   show_percentage = FALSE,
#   fill_color = mycolor[1:length(x)],
#   #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4
# )
# 
# ggsave(plot = venn_plot, filename = 'runtime/8.venn/venn_all.pdf', width = 6, height = 6)
# 
# 
# library(VennDiagram)
# 
# v <- venn.diagram(x[1:6], filename = NULL)
# grid.newpage()
# grid.draw(v)

load('runtime/4.DEM.new/res_DEM.rda')

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}
contrasts_target.treatment <- c('eHx-pHx:18', 'eHx-pHx:36', 'eHx-pHx:72')

vip1 <- readRDS('runtime/6.VIP_score/Treatment.18/vip.mat.rds')
vip2 <- readRDS('runtime/6.VIP_score/Treatment.36/vip.mat.rds')
vip3 <- readRDS('runtime/6.VIP_score/Treatment.72/vip.mat.rds')

compounds.18 <- res_DEM$Index[which(res_DEM$contrast %in% 'eHx-pHx:18' & 
                                      res_DEM$Index %in% rownames(vip1)[which(rowMeans(vip1[, 1:2]) > 1.5)] &
                                      res_DEM$q < 0.05)]
compounds.36 <- res_DEM$Index[which(res_DEM$contrast %in% 'eHx-pHx:36' & 
                                      res_DEM$Index %in% rownames(vip1)[which(rowMeans(vip1[, 1:2]) > 1.5)] &
                                      res_DEM$q < 0.05)]
compounds.72 <- res_DEM$Index[which(res_DEM$contrast %in% 'eHx-pHx:72' & 
                                      res_DEM$Index %in% rownames(vip1)[which(rowMeans(vip1[, 1:2]) > 1.5)] &
                                      res_DEM$q < 0.05)]



df_plot <- data.frame(
  compound = c(compounds.18, compounds.36, compounds.72),
  contrast = c(rep('eHx-pHx:18', length(compounds.18)),
               rep('eHx-pHx:36', length(compounds.36)),
               rep('eHx-pHx:72', length(compounds.72))),
  stringsAsFactors = F
)

x <- split_tibble(df_plot, 'contrast')

mycolor <- c("blue", "yellow", "steelblue", "red", 'salmon', "green")

venn_plot <- ggvenn(
  x, 
  show_percentage = FALSE,
  fill_color = mycolor[1:length(x)],
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave(plot = venn_plot, filename = 'runtime/8.venn/venn.compared_treatment.in.3time_point.pdf', width = 6, height = 4.5)






