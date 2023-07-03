
load('runtime/2.PCA/mSet.Rdata')

pc1 = mSet$analSet$pca$x[, 1]
pc2 = mSet$analSet$pca$x[, 2]

xlabel = paste0("PC ", 1, " (", round(100 * mSet$analSet$pca$variance[1], 
                                     1), "%)")
ylabel = paste0("PC ", 2, " (", round(100 * mSet$analSet$pca$variance[2], 
                                     1), "%)")


df_plot <- data.frame(X = pc1, Y = pc2)

df_plot$Group <- gsub('^([^-]+)-.*$', '\\1', rownames(df_plot), perl = T)
df_plot$Time <- gsub('^.*-(\\d+)-.*$', '\\1', rownames(df_plot), perl = T)
df_plot$Time[grep('sham', df_plot$Time)] <- '0'
df_plot$sample <- rownames(df_plot)

pca_2d <- ggplot(df_plot, aes(X, Y)) +
  stat_ellipse(geom = "polygon", type = 'norm', alpha = 0.5, aes(fill = Time)) +
  geom_point(aes(shape=Group), size = 2) +
  geom_text(aes(label=sample), size = 3, nudge_y = 1) +
  scale_shape_manual(values = c('eHx' = 0, 'pHx' = 1, 'sham' = 2)) +
  labs(x = xlabel,
       y = ylabel,
       title = 'Scores Plot') +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = .5, size = 15)
  )

ggsave(plot = pca_2d, filename = 'runtime/2.PCA/pca_score2d_0_dpi72.new.pdf', width = 7.5, height = 6.5)



