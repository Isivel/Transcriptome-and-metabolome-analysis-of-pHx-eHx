
library(ggvenn)

outDir <- './runtime/3.venn'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)


load('runtime/2.DEG/res_DEGs.condition.rda')

contrasts <- rbind(c(0, 1, 0, 0, 0, 1, 0, 0), # eHx vs pH in 18
                   c(0, 1, 0, 0, 0, 0, 1, 0), # eHx vs pH in 36
                   c(0, 1, 0, 0, 0, 0, 0, 1) # eHx vs pH in 72
)

res_DEGs.condition$contrast[res_DEGs.condition$contrast == '01000100'] <- 'eHx-pHx:18'
res_DEGs.condition$contrast[res_DEGs.condition$contrast == '01000010'] <- 'eHx-pHx:36'
res_DEGs.condition$contrast[res_DEGs.condition$contrast == '01000001'] <- 'eHx-pHx:72'



df_plot <- res_DEGs.condition[res_DEGs.condition$padj < 0.001 & res_DEGs.condition$log2FoldChange > 1, c('geneID', 'contrast')]

split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}
x <- split_tibble(df_plot, 'contrast')

mycolor <- c("blue", "yellow", "steelblue", "red", 'salmon', "green")

venn_plot <- ggvenn(
  x, 
  show_percentage = FALSE,
  fill_color = mycolor[1:length(x)],
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave(plot = venn_plot, filename = 'runtime/3.venn/venn_all.pdf', width = 6, height = 4.5)



