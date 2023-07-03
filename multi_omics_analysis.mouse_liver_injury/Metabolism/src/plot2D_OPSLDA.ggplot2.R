
load('runtime/3.OPLSDA/mSet.Rdata')

lv1 <- mSet$analSet$oplsda$scoreMN[, 1]
lv2 <- mSet$analSet$oplsda$orthoScoreMN[, 1]
xlabel <- paste("T score [1]", "(", round(100 * mSet$analSet$oplsda$modelDF["p1", 
                                                                               "R2X"], 1), "%)")
ylabel <- paste("Orthogonal T score [1]", "(", round(100 * 
                                                       mSet$analSet$oplsda$modelDF["o1", "R2X"], 1), "%)")
text.lbls <- substr(rownames(mSet$df_plotSet$norm), 1, 12)
lvs <- levels(mSet$df_plotSet$cls)


df_plot <- data.frame(X = lv1, Y = lv2)

df_plot$Group <- gsub('^([^-]+)-.*$', '\\1', rownames(df_plot), perl = T)
df_plot$Time <- gsub('^.*-(\\d+)-.*$', '\\1', rownames(df_plot), perl = T)
df_plot$Time[grep('sham', df_plot$Time)] <- '0'

opls_2d <- ggplot(df_plot, aes(X, Y)) +
  stat_ellipse(geom = "polygon", alpha = 0.5,aes(fill = Time)) +
  geom_point(aes(shape=Group), size = 2) +
  #geom_text(aes(label=X), size = 3, nudge_y = 1) +
  scale_shape_manual(values = c('eHx' = 0, 'pHx' = 1, 'sham' = 2)) +
  labs(x = xlabel,
       y = ylabel,
       title = 'Scores Plot') +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = .5, size = 15)
  )

ggsave(plot = opls_2d, filename = 'runtime/3.OPLSDA/opls_score2d_0_dpi72.new.pdf', width = 7.5, height = 6.5)



