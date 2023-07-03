
load('runtime/3.OPLSDA/mSet.Rdata')

mSetObj <- mSet
imgName = 'runtime/3.OPLSDA/opls_score2d_0_dpi72.new.pdf'
width <- NA
if (is.na(width)) {
  w <- 9
}else if (width == 0) {
  w <- 7.2
}else {
  w <- width
}
h <- w

Cairo::Cairo(file = imgName, unit = "in", dpi = 72, width = w, 
             height = h, type = 'pdf', bg = "white")
par(mar = c(5, 5, 3, 3))
lv1 <- mSetObj$analSet$oplsda$scoreMN[, 1]
lv2 <- mSetObj$analSet$oplsda$orthoScoreMN[, 1]
xlabel <- paste("T score [1]", "(", round(100 * mSetObj$analSet$oplsda$modelDF["p1", 
                                                                               "R2X"], 1), "%)")
ylabel <- paste("Orthogonal T score [1]", "(", round(100 * 
                                                       mSetObj$analSet$oplsda$modelDF["o1", "R2X"], 1), "%)")
text.lbls <- substr(rownames(mSetObj$dataSet$norm), 1, 12)
lvs <- levels(mSetObj$dataSet$cls)
pts.array <- array(0, dim = c(100, 2, length(lvs)))
for (i in 1:length(lvs)) {
  inx <- mSetObj$dataSet$cls == lvs[i]
  groupVar <- var(cbind(lv1[inx], lv2[inx]), na.rm = T)
  groupMean <- cbind(mean(lv1[inx], na.rm = T), mean(lv2[inx], 
                                                     na.rm = T))
  pts.array[, , i] <- ellipse::ellipse(groupVar, centre = groupMean, 
                                       level = reg, npoints = 100)
}
xrg <- range(lv1, pts.array[, 1, ])
yrg <- range(lv2, pts.array[, 2, ])
x.ext <- (xrg[2] - xrg[1])/12
y.ext <- (yrg[2] - yrg[1])/12
xlims <- c(xrg[1] - x.ext, xrg[2] + x.ext)
ylims <- c(yrg[1] - y.ext, yrg[2] + y.ext)
cols <- GetColorSchema(mSetObj$dataSet$cls, grey.scale == 
                         1)
uniq.cols <- unique(cols)
plot(lv1, lv2, xlab = xlabel, xlim = xlims, ylim = ylims, 
     ylab = ylabel, type = "n", main = "Scores Plot")
grid(col = "lightgray", lty = "dotted", lwd = 1)
legend.nm <- unique(as.character(mSetObj$dataSet$cls))
if (length(uniq.cols) > 1) {
  names(uniq.cols) <- legend.nm
}
for (i in 1:length(lvs)) {
  if (length(uniq.cols) > 1) {
    polygon(pts.array[, , i], col = adjustcolor(uniq.cols[lvs[i]], 
                                                alpha = 0.2), border = NA)
  }
  else {
    polygon(pts.array[, , i], col = adjustcolor(uniq.cols, 
                                                alpha = 0.2), border = NA)
  }
  if (grey.scale) {
    lines(pts.array[, , i], col = adjustcolor("black", 
                                              alpha = 0.5), lty = 2)
  }
}
pchs <- GetShapeSchema(mSetObj, show, grey.scale)
if (grey.scale) {
  cols <- rep("black", length(cols))
}
if (show == 1) {
  text(lv1, lv2, label = text.lbls, pos = 4, xpd = T, cex = 0.75)
  points(lv1, lv2, pch = pchs, col = cols)
}else {
  if (length(uniq.cols) == 1) {
    points(lv1, lv2, pch = pchs, col = cols, cex = 1)
  }
  else {
    if (grey.scale == 1 | (exists("shapeVec") && all(shapeVec >= 
                                                     0))) {
      my.cols <- adjustcolor(cols, alpha.f = 0.4)
      my.cols[pchs == 21] <- "black"
      points(lv1, lv2, pch = pchs, col = my.cols, bg = adjustcolor(cols, 
                                                                   alpha.f = 0.4), cex = 1.8)
    }
    else {
      points(lv1, lv2, pch = 21, bg = adjustcolor(cols, 
                                                  alpha.f = 0.4), cex = 2)
    }
  }
}
uniq.pchs <- unique(pchs)
if (grey.scale) {
  uniq.cols <- "black"
}
legend("topright", legend = legend.nm, pch = uniq.pchs, col = uniq.cols)
dev.off()
return(.set.mSet(mSetObj))




