
outDir <- './runtime/5.S_plot'
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

load('runtime/3.OPLSDA/mSet.Rdata')

res_sig <- readRDS('runtime/4.DEM/res_2wayAnova.rds')

splot_dat <- as.data.frame(mSet$analSet$oplsda$splot.mat)
vip <- mSet$analSet$oplsda$vipVn

splot_dat$sig <- 'N'
splot_dat$sig[which(rownames(splot_dat) %in% names(vip)[vip > 1] & 
                      rownames(splot_dat) %in% res_sig$me[which(res_sig$q.time < 0.05)])] <- 'Y'

colnames(splot_dat)[2:3] <- c('p1', 'pcoor1') 

s_plot <- ggplot(data = splot_dat, aes(x = p1, y = pcoor1)) + #'p[1]', y = 'p(corr)[1]'
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c('N' = 'grey', 'Y' = 'salmon')) +
  labs(x = 'p[1]', y = 'p(corr)[1]') +
  theme_bw() +
  theme(
    legend.position = 'none'
  )
ggsave(plot = s_plot, width = 7.5, height = 5, filename = 'runtime/5.S_plot/S_plot.pdf')


# #### boxplot
# dat.met.norm <- readRDS('runtime/1.data/all_sample.norm.RDS')
# dat.sig <- dat.met.norm[which(rownames(dat.met.norm) %in% names(vip)[order(vip, decreasing = T)[1:30]] & 
#                                 rownames(dat.met.norm) %in% res_sig$me[which(res_sig$q.time < 0.05)]), ]
# dat.sig <- as.data.frame(t(dat.sig))
# dat.sig$sample <- rownames(dat.sig)
# dat.sig_df <- melt(dat.sig)
# 
# dat.sig_df$Group <- gsub('^([^-]+)\\-.*$', '\\1', dat.sig_df$sample, perl = T)
# dat.sig_df$Time <- gsub('^.*\\-(\\d+)\\-.*$', '\\1', dat.sig_df$sample, perl = T)
# dat.sig_df$Time[grep('sham', dat.sig_df$Time)] <- '0'
# 
# ggplot(data = dat.sig_df, aes(x = Time, y = value)) +
#   geom_boxplot(aes(fill = Group)) +
#   theme_bw() +
#   facet_wrap(~variable, ncol = 5)


#### change coordinate axis
# theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
#                            color = "black", size = 1, 
#                            xlab = "x", ylab = "y",
#                            ticks = 5,
#                            textsize = 3,
#                            xlimit = max(abs(xvals)),
#                            ylimit = max(abs(yvals)),
#                            epsilon = max(xlimit,ylimit)/50){
#   
#   #INPUT:
#   #xvals .- Values of x that will be plotted
#   #yvals .- Values of y that will be plotted
#   #xgeo  .- x intercept value for y axis
#   #ygeo  .- y intercept value for x axis
#   #color .- Default color for axis
#   #size  .- Line size for axis
#   #xlab  .- Label for x axis
#   #ylab  .- Label for y axis
#   #ticks .- Number of ticks to add to plot in each axis
#   #textsize .- Size of text for ticks
#   #xlimit .- Limit value for x axis 
#   #ylimit .- Limit value for y axis
#   #epsilon .- Parameter for small space
#   
#   
#   #Create axis 
#   xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
#   yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))
#   
#   #Add axis
#   theme.list <- 
#     list(
#       theme_void(), #Empty the current theme
#       geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
#       geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
#       annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
#       annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
#       xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
#       ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
#     )
#   
#   #Add ticks programatically
#   ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
#   ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)
#   
#   #Add ticks of x axis
#   nlist <- length(theme.list)
#   for (k in 1:ticks){
#     
#     #Create data frame for ticks in x axis
#     xtick <- data.frame(xt = rep(ticks_x[k], 2), 
#                         yt = c(xgeo + epsilon, xgeo - epsilon))
#     
#     #Create data frame for ticks in y axis
#     ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
#                         yt = rep(ticks_y[k], 2))
#     
#     #Add ticks to geom line for x axis
#     theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
#                                              data = xtick, size = size, 
#                                              color = color)
#     
#     #Add labels to the x-ticks
#     theme.list[[nlist + 4*k-2]] <- annotate("text", 
#                                             x = ticks_x[k], 
#                                             y = ygeo - 2.5*epsilon,
#                                             size = textsize,
#                                             label = paste(ticks_x[k]))
#     
#     
#     #Add ticks to geom line for y axis
#     theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
#                                              data = ytick, size = size, 
#                                              color = color)
#     
#     #Add labels to the y-ticks
#     theme.list[[nlist + 4*k]] <- annotate("text", 
#                                           x = xgeo - 2.5*epsilon, 
#                                           y = ticks_y[k],
#                                           size = textsize,
#                                           label = paste(ticks_y[k]))
#   }
#   
#   #Add theme
#   #theme.list[[3]] <- 
#   return(theme.list)
# }
# 
# ggplot(splot_dat, aes(x = p1, y = pcoor1)) +
#   theme_geometry(splot_dat$p1, splot_dat$pcoor1) +
#   geom_point(aes(color = sig)) +
#   scale_color_manual(values = c('N' = 'grey', 'Y' = 'salmon')) +
#   labs(x = 'p[1]', y = 'p(corr)[1]')
