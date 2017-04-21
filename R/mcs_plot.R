McsPlot <- function(lsmeans, pval, trmts, Edescriptors, Emethods,
                    single.desc, metric) {

  # Re-order for plotting
  # large lsmeans are good, unless using RMSE
  if (!metric %in% c("RMSE", "error rate")) {
    o <- order(lsmeans, decreasing = TRUE)
  } else {
    o <- order(lsmeans, decreasing = FALSE)
  }
  lsm <- lsmeans[o]
  pval <- pval[o, o]
  trmts <- trmts[o]

  # Emethods <-
  # c('ENet','RF','KNN','LAR','NNet','PCR','PLS','PLSLDA','RPart','Ridge','SVM','Tree')
  n <- length(lsm)

  # TO DO: Is there a better way to do change figure margins?
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(pty = "s", oma = c(0, 0, 0, 0), mar = c(5, 5.5, 6, 5.5) + 0,
      mfrow = c(1,1), mgp = c(3, 1, 0))

  plot(1:n, 1:n, ylim = c(1, n), xaxt = "n", yaxt = "n", pty = "s", col = 0, xlab = "",
       ylab = "")
  # axis( side=1, at=1:n, labels=trmts, cex.axis=.5, las=3, tck=1, lty=0,
  # col='grey', line=-.6)
  axis(side = 1, at = 1:n, labels = round(lsm, digits = 2), cex.axis = 0.5, las = 3,
       tck = 1, lty = 0, col = "grey", line = -0.6)
  if (single.desc)
    axis(side = 1, at = 1:n, labels = Emethods[trmts%%100], cex.axis = 0.5, las = 3,
         tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5) else {
           axis(side = 1, at = 1:n, labels = Emethods[trmts%%100], cex.axis = 0.5, las = 3,
                tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
           axis(side = 1, at = 1:n, labels = Edescriptors[trmts%/%100], cex.axis = 0.5,
                las = 3, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 2.3)
         }
  # axis( side=2, at=1:n, labels=trmts[n:1], cex.axis=.5, las=1, tck=1, lty=0,
  # col='grey', line=-.6)
  axis(side = 2, at = 1:n, labels = round(lsm[n:1], digits = 2), cex.axis = 0.5,
       las = 1, tck = 1, lty = 0, col = "grey", line = -0.6)
  if (single.desc)
    axis(side = 2, at = 1:n, labels = Emethods[trmts[n:1]%%100], cex.axis = 0.5,
         las = 1, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
  else {
    axis(side = 2, at = 1:n, labels = Emethods[trmts[n:1]%%100], cex.axis = 0.5,
         las = 1, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
    axis(side = 2, at = 1:n, labels = Edescriptors[trmts[n:1]%/%100], cex.axis = 0.5,
         las = 1, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 2.3)
  }
  axis(side = 3, at = 1:n, labels = round(lsm, digits = 2), cex.axis = 0.5, las = 3,
       tck = 1, lty = 0, col = "grey", line = -0.6)
  if (single.desc)
    axis(side = 3, at = 1:n, labels = Emethods[trmts%%100], cex.axis = 0.5, las = 3,
         tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
  else {
    axis(side = 3, at = 1:n, labels = Edescriptors[trmts%/%100], cex.axis = 0.5,
         las = 3, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
    axis(side = 3, at = 1:n, labels = Emethods[trmts%%100], cex.axis = 0.5, las = 3,
         tck = 1, lty = 0, col = "grey", tick = FALSE, line = 1.7)
  }
  axis(side = 4, at = 1:n, labels = round(lsm[n:1], digits = 2), cex.axis = 0.5,
       las = 1, tck = 1, lty = 0, col = "grey", line = -0.6)
  if (single.desc)
    axis(side = 4, at = 1:n, labels = Emethods[trmts[n:1]%%100], cex.axis = 0.5,
         las = 1, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
  else {
    axis(side = 4, at = 1:n, labels = Edescriptors[trmts[n:1]%/%100], cex.axis = 0.5,
         las = 1, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 0.5)
    axis(side = 4, at = 1:n, labels = Emethods[trmts[n:1]%%100], cex.axis = 0.5,
         las = 1, tck = 1, lty = 0, col = "grey", tick = FALSE, line = 1.7)
  }
  pval[is.na(pval)] <- 1
  plot.pval <- pval
  for (i in 1:n) {
    for (j in 1:n) {
      if (pval[i, j] >= 0.9)
        plot.pval[i, j] <- 4
      else if (pval[i, j] >= 0.1)
        plot.pval[i, j] <- 3
      else if (pval[i, j] >= 0.05)
        plot.pval[i, j] <- 2
      else if (pval[i, j] >= 0.01)
        plot.pval[i, j] <- 1
      else
        plot.pval[i, j] <- 5
    }
  }
  image(x = 1:n, y = 1:n, z = plot.pval[1:n, n:1], col = c("#ffffcc", "#a1dab4",
                                                           "#41b6c4", "#253494", "white")[sort(unique(c(plot.pval)))], add = TRUE)
  abline(h = 0.5:(n + 0.5), col = "grey")
  abline(v = 0.5:(n + 0.5), col = "grey")
  mtext(side = 3, line = 5, text = "Multiple Comparisons Similarity (MCS) Plot",
        cex = 0.8)
  mtext(side = 1, line = -4.75, adj = -0.04, cex = 0.65,
        text = "                Model\n               Performance-> ",
        outer = TRUE)
  if (single.desc)
    mtext(side = 1, line = -3.5, adj = -0.04, cex = 0.65,
          text = "                Method-> ", outer = TRUE)
  else mtext(side = 1, line = -2.1, adj = -0.04, cex = 0.65,
             text = "                Method-> \n \n                Descriptor->",
             outer = TRUE)
  mtext(side = 3, line = -5.45, adj = 1.04, cex = 0.65,
        text = "Model                 ",
        outer = TRUE)
  mtext(side = 3, line = -5.9, adj = 1.04, cex = 0.65,
        text = " <-Performance                ",
        outer = TRUE)
  if (single.desc)
    mtext(side = 3, line = -4.6, adj = 1.04, cex = 0.65,
          text = " <-Method                ",
          outer = TRUE)
  else mtext(side = 3, line = -4.6, adj = 1.04, cex = 0.65,
             text = " <-Method                \n \n <-Descriptor                ",
             outer = TRUE)
  temp <- legend("bottomleft", legend = c(expression(paste("         p-val<.01;  ",
                                                           0.01 <= alpha)), expression(paste(0.01 <= p, "-val<.05;  ", 0.05 <= alpha)),
                                          expression(paste(0.05 <= p, "-val<.10;    ", 0.1 <= alpha)), expression(paste("  ",
                                                                                                                        0.1 <= p, "-val<.90;    ", 0.9 <= alpha, "      ")), expression(paste("  ",
                                                                                                                                                                                              0.9 <= p, "-val<1"))), title = "\n", fill = c("white", "#ffffcc", "#a1dab4",
                                                                                                                                                                                                                                            "#2c7fb8", "#253494"), bg = "white", y.intersp = 0.9, cex = 0.55, inset = -0.01)
  temp2 <- legend(temp$rect$left - 0.35, temp$rect$top, legend = c("Multiplicity-adjusted p-values;",
                                                                   expression(paste("QSAR models different at level ", alpha, "  ", sep = ""))),
                  cex = 0.52, bg = "white", y.intersp = 0.9, plot = F)
  legend(temp$rect$left, temp$rect$top + temp2$rect$h, legend = c("Multiplicity-adjusted p-values;",
                                                                  expression(paste("QSAR models different at level ", alpha, "  ", sep = ""))),
         cex = 0.52, bg = "white", y.intersp = 0.9)
  box()
}
