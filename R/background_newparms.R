
#--------------------------------------------------------------------------------

HitMapPages <- function(max.per.page, x, c, hc, title, ...) {
  if (attr(hc, "members") <= max.per.page)
    heatmap(as.matrix(x[order.dendrogram(hc), ]), col = c, Rowv = hc, Colv = NA,
      scale = "none", main = title, ...)
  else if (attr(hc, "height") == 0) {
    done <- 0
    while (!(done)) {
      if (nrow(x) > max.per.page) {
        mx <- ceiling(nrow(x)/max.per.page)
        mxdone <- 0
        max.per.page <- max.per.page - 1
        while (!(mxdone)) {
          if (ceiling(nrow(x)/max.per.page) > mx) {
          mxdone <- 1
          max.per.page <- max.per.page + 1
          } else max.per.page <- max.per.page - 1
        }
        heatmap(as.matrix(x[(1:max.per.page), ]), col = c, Rowv = NA, Colv = NA,
          scale = "none", main = title, ...)
        x <- x[-(1:max.per.page), ]
      } else {
        heatmap(as.matrix(x), col = c, Rowv = NA, Colv = NA, scale = "none",
          main = title, ...)
        done <- 1
      }
    }
  } else {
    hc.cut <- cut(hc, 2)
    HitMapPages(max.per.page, x, c, hc.cut$lower[[1]], title, ...)
    HitMapPages(max.per.page, x, c, hc.cut$lower[[2]], title, ...)
  }
}

HitMap <- function(yhat.list, y, max.select = NA, max.per.page = 1000, max.cont.y = 50,
  labels = NA, y.labels = NA, title = "", descriptors, ...) {

  if (missing(labels))
    labels <- names(yhat.list)
  if (missing(y.labels))
    y.labels <- row.names(yhat.list)

  if (sum(!(y %in% c(1, 0))))
    classify <- "N" else classify <- "Y"

  # if max.select not specified use min(300,#compounds/4)
  if (missing(max.select))
    max.select <- min(300, (length(y)/4))

  if (sum(!(y %in% c(1, 0))))
    classify <- "N" else classify <- "Y"
  list.index <- 1
  while (list.index <= length(labels)) {
    if (length(labels) == 1)
      yhat <- unlist(yhat.list) else yhat <- unlist(yhat.list[list.index])
    yhat <- rank(-yhat)
    if (classify == "Y")
      yhat <- yhat[y == 1] else yhat <- yhat[rank(-y, ties = "rand") <= max.cont.y]
    if (list.index == 1)
      x <- data.frame(yhat) else x <- data.frame(x, yhat)
    list.index <- list.index + 1
  }
  if (length(labels) == 1) {
    labels <- c(labels, "")
    x <- data.frame(x, rep(max.select, length(yhat)))
  }
  names(x) <- labels
  if (classify == "Y")
    row.names(x) <- y.labels[y == 1]
  else row.names(x) <- y.labels[rank(-y, ties = "rand") <= max.cont.y]

  for (i in 1:nrow(x)) for (j in 1:ncol(x)) if (x[i, j] > max.select)
    x[i, j] <- max.select

  c <- c()
  for (i in 1:max.select) c <- c(c, rgb((max.select - (0.1 * i)), (0.9 * i), (0.9 * i),
                                        maxColorValue = max.select))

  if ((missing(descriptors)) || (nrow(descriptors) == 0)) {
    x <- x[order(rowSums(x != max.select)), order(colSums(x != max.select))]
    heatmap(as.matrix(x), col = c, Rowv = NA, Colv = NA,
            scale = "none", main = title, ...)
  } else {
    if (classify == "Y")
      descriptors <- descriptors[y == 1, ]
    else descriptors <- descriptors[rank(-y, ties = "rand") <= max.cont.y, ]
    x <- x[, order(colSums(x != max.select))]
    d <- dist(descriptors, "bin")
    hc <- as.dendrogram(hclust(d))
    if (attr(hc, "members") <= max.per.page)
      heatmap(as.matrix(x), col = c, Rowv = hc, Colv = NA, scale = "none",
        main = title, ...) else HitMapPages(max.per.page, x, c, hc, title, ...)
  }

}

plotSingleRun <- function(cml.result, splitidx, desidx) {
  y <- cml.result$response
  desc <- cml.result$data[[desidx]]

  # is this necessary, wasnt data cleaning already done?
  y <- y[apply(desc, 1, function(x) sum(is.na(x)) == 0)]
  desc <- subset(desc, apply(desc, 1, function(x) sum(is.na(x)) == 0))
  desc <- subset(desc, !is.na(y))
  y <- y[!is.na(y)]
  if (sum(!(y %in% c(1, 0))))
    classify <- "N" else classify <- "Y"
  preds <- cml.result$all.preds[[splitidx]][[desidx]]
  if (classify == "Y")
    probs <- cml.result$all.probs[[splitidx]][[desidx]]
  # bitmap(paste(plotfile,'.bmp',sep=''))
  # pdf(paste(plotfile,'.pdf',sep=''))
  if ((classify == "Y") && (ncol(probs) > 1)) {
    HitCurve(probs[, -1], y = probs[, 1], phat.labels = names(probs)[-1])
    ContCurve(preds[, !(names(preds) %in% names(probs))],
              yhat.labels = names(preds)[!(names(preds) %in%
                                             names(probs))], y = preds[, 1],
              curves.only = TRUE, start.col = (ncol(probs) - 1))
  } else {
    ContCurve(preds[, -1], y = preds[, 1], yhat.labels = names(preds)[-1])
  }
  # dev.off()
  # pdf(paste(plotfile,'s.pdf',sep=''))
  if ((classify == "Y") && (ncol(probs) > 1)) {
    HitCurve(probs[, -1], y = probs[, 1], phat.labels = names(probs)[-1])
    ContCurve(preds[, !(names(preds) %in% names(probs))],
              yhat.labels = names(preds)[!(names(preds) %in%
                                             names(probs))], y = preds[, 1],
              curves.only = TRUE, start.col = (ncol(probs) - 1))
    # HitMap(data.frame(probs[,-1],preds[,!(names(preds) %in% names(probs))]),
    # y=preds[,1], y.labels=row.names(desc), descriptors=desc,
    # labels=c(names(probs)[-1],names(preds)[!(names(preds) %in% names(probs))]))
  } else {
    ContCurve(preds[, -1], y = preds[, 1], yhat.labels = names(preds)[-1])
    # HitMap(preds[,-1], y=preds[,1], y.labels=row.names(desc), descriptors=desc,
    # labels=names(preds)[-1])
  }
  # dev.off()
}


#--------------------------------------------------------------------------------
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

  par(pty = "s", oma = c(0, 0, 0, 0), mar = c(5, 5.5, 6, 5.5) + 0)
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
        text = "                 Model\n                Accuracy-> ",
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
        text = " <-Accuracy               ",
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

#--------------------------------------------------------------------------------
# This section contains code for combining multiple descriptor sets
#--------------------------------------------------------------------------------

# combinedescriptors <- function(bndir,phdir,apdir,fpdir,chdir,response,dest,root,num.actives=250) {
#   tryCatch( combinedescriptors.try(bndir,phdir,apdir,fpdir,chdir,response,dest,root,num.actives=250),
#             error = function(e) {warning(paste("ERROR...Descriptors not combined:",e$message))} )
# }

plot.ChemModLab <- function(cml.result) {
  # bndir,phdir,apdir,fpdir,chdir: contain the dir name of the various descriptor
  #                                set runs, e.g. Split1/BurdenNumbers/2007-02-26-16-18-07,
  #                                note no begining or ending slashes, any of these
  #                                may be missing, at least one required
  # response: the name of the response file, e.g Outcome.csv
  # dest: the dest directory, where the combined files will reside, no
  #       begining or ending slashes, e.g Split1
  # root: the root dir, ends in a slash
  # num.actives: used for continuous responses where the number of actives
  #              can't be easily determined

#   if (missing("bndir")) bndir=""
#   if (missing("phdir")) phdir=""
#   if (missing("apdir")) apdir=""
#   if (missing("fpdir")) fpdir=""
#   if (missing("chdir")) chdir=""
#
#   i<-1
#   preds<-list()
#   if (bndir!="") { preds[[i]]<-read.csv(paste(root,bndir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
#   i<-i+1 }
#   if (phdir!="") { preds[[i]]<-read.csv(paste(root,phdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
#   i<-i+1 }
#   if (apdir!="") { preds[[i]]<-read.csv(paste(root,apdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
#   i<-i+1 }
#   if (fpdir!="") { preds[[i]]<-read.csv(paste(root,fpdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
#   i<-i+1 }
#   if (chdir!="") { preds[[i]]<-read.csv(paste(root,chdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
#   i<-i+1 }

  nsplit <- length(cml.result$all.preds)
  for (splidx in 1:nsplit) {
    preds <- cml.result$all.preds[[splidx]]
    titles <- paste0("Split ", splidx, " : ", gsub("_", " ", names(preds)))
    abrev.titles <- c()
    for (j in 1:length(preds)) {
      abrev.titles <- c(abrev.titles, paste0("Des", j))
    }
  #   if (bndir!="") { titles<-c(titles,"Burden Numbers")
  #   abrev.titles<-c(abrev.titles,"B#") }
  #   if (phdir!="") { titles<-c(titles,"Pharmacophores")
  #   abrev.titles<-c(abrev.titles,"Ph") }
  #   if (apdir!="") { titles<-c(titles,"Atom Pairs")
  #   abrev.titles<-c(abrev.titles,"AP") }
  #   if (fpdir!="") { titles<-c(titles,"Fragment Pairs")
  #   abrev.titles<-c(abrev.titles,"FP") }
  #   if (chdir!="") { titles<-c(titles,"Carhart Atom Pairs")
  #   abrev.titles<-c(abrev.titles,"CAP") }

    y <- cml.result$responses
    if (sum(!(y %in% c(1, 0))))
      classify <- "N" else classify <- "Y"
    if (classify == "Y")
      num.actives <- sum(y)
    if (classify == "Y")
      probs <- cml.result$all.probs[[splidx]]

    # TO DO why removing KNN, etc predictions?
    if (classify == "N") {
      for (i in length(titles))
        preds[[i]] <- preds[[i]][, !(names(preds[[i]]) %in% c("KNN", "NNet", "PLSLDA"))]
    }

    # desc<-data.frame()
#     if (chdir!="")
#       desc<-read.csv(paste(root,"Carharts.csv",sep=""),header=TRUE, row.names=1)
#
    # pdf(paste(root,dest,"/plots.pdf",sep=""))
    #Plot all methods for each descriptor set
    for (i in 1:length(titles)) {
      if ((classify == "Y") && (ncol(probs[[i]]) > 1)) {
        HitCurve(probs[[i]][, -1], y = y, title = titles[i], phat.labels = names(probs[[i]])[-1])
        ContCurve(preds[[i]][, !(names(preds[[i]]) %in% names(probs[[i]]))],
          y = y, curves.only = TRUE, start.col = (ncol(probs[[i]]) - 1),
          title = titles[i], yhat.labels = names(preds[[i]])[!(names(preds[[i]]) %in%
          names(probs[[i]]))])
#         if (ncol(desc)==0)
#           hit.map(data.frame(probs[[i]][,-1],preds[[i]][,!(names(preds[[i]]) %in% names(probs[[i]]))]),
#                   y=y, title=titles[i], labels=c(names(probs[[i]])[-1],names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))]))
#         else
#           hit.map(data.frame(probs[[i]][,-1],preds[[i]][,!(names(preds[[i]]) %in% names(probs[[i]]))]),
#                   y=y, y.labels=row.names(desc), descriptors=desc, title=titles[i],
#                   labels=c(names(probs[[i]])[-1],names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))]))

      } else {
        ContCurve(preds[[i]][, -1], y = y, title = titles[i], yhat.labels = names(preds[[i]])[-1])
#         if (ncol(desc)==0)
#           hit.map(preds[[i]][,-1], y=y, labels=names(preds[[i]])[-1],
#                   title=titles[i], max.cont.y=num.actives, max.per.page=num.actives)
#         else
#           hit.map(preds[[i]][,-1], y=y, y.labels=row.names(desc), descriptors=desc,
#                   title=titles[i], max.cont.y=num.actives, max.per.page=num.actives, labels=names(preds[[i]])[-1])

      }
    }

    # Plot all descriptors for each method
    meths <- c()
    for (i in length(titles)) meths <- c(meths, names(preds[[i]])[-1])
    meths <- unique(meths)
    if (classify == "Y") {
      pmeths <- c()
      for (i in length(titles)) pmeths <- c(pmeths, names(probs[[i]])[-1])
      pmeths <- unique(pmeths)
      meths <- meths[!(meths %in% pmeths)]
      for (j in 1:length(pmeths)) {
        p <- data.frame()
        p.labels <- c()
        for (i in 1:length(titles)) {
          if (sum(names(probs[[i]]) %in% pmeths[j]) == 1) {
          if (ncol(p) == 0)
            p <- data.frame(probs[[i]][, pmeths[j]])
          else p <- data.frame(p, probs[[i]][, pmeths[j]])
          names(p)[ncol(p)] <- abrev.titles[i]
          p.labels <- c(p.labels, abrev.titles[i])
          }
        }
#         if (ncol(p)>0) {
#           hit.curve(p, y=y, title=pmeths[j], phat.labels=p.labels)
#           if (ncol(desc)==0)
#             hit.map(p, y=y, title=pmeths[j], cexCol=1, labels=p.labels)
#           else
#             hit.map(p, y=y, y.labels=row.names(desc), descriptors=desc, title=pmeths[j], cexCol=1, labels=p.labels)
#         }
      }
    }
    for (j in 1:length(meths)) {
      p <- data.frame()
      p.labels <- c()
      for (i in 1:length(titles)) {
        if (sum(names(preds[[i]]) %in% meths[j]) == 1) {
          if (ncol(p) == 0)
          p <- data.frame(preds[[i]][, meths[j]]) else p <- data.frame(p, preds[[i]][, meths[j]])
          names(p)[ncol(p)] <- abrev.titles[i]
          p.labels <- c(p.labels, abrev.titles[i])
        }
      }
      if (ncol(p) > 0) {
        ContCurve(p, y = y, title = paste("Split", splidx, ":", meths[j]),
          yhat.labels = p.labels)
#         if (ncol(desc)==0)
#           hit.map(p, y=y, title=meths[j],
#                   cexCol=1, max.cont.y=num.actives, max.per.page=num.actives, labels=p.labels)
#         else
#           hit.map(p, y=y, y.labels=row.names(desc), descriptors=desc, title=meths[j],
#                   cexCol=1, max.cont.y=num.actives, max.per.page=num.actives, labels=p.labels)
      }
    }

#     #Create plot.pdf
#     i<-length(titles)
#     if ((classify=="Y") && (ncol(probs[[i]])>1)) {
#       hit.curve(probs[[i]][,-1], y=y, title=titles[i], phat.labels=names(probs[[i]])[-1])
#       cont.curve(preds[[i]][,!(names(preds[[i]]) %in% names(probs[[i]]))],
#                  y=y, curves.only=TRUE, start.col=(ncol(probs[[i]])-1), title=titles[i],
#                  yhat.labels=names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))])
#     }
#     else {
#       cont.curve(preds[[i]][,-1], y=y, title=titles[i], yhat.labels=names(preds[[i]])[-1])
#     }
  }
    # dev.off()

#     #Create combined predictions
#     head<-","
#     p<-data.frame()
#     for (i in 1:length(titles)) {
#       head<-paste(head,",",titles[i],sep="")
#       if (ncol(preds[[i]])>2)
#         for (j in 2:(ncol(preds[[i]])-1))
#           head<-paste(head,",",sep="")
#         if (ncol(p)==0) p<-data.frame(preds[[i]])
#         else p<-data.frame(p,subset(preds[[i]],select=-1),check.names=FALSE)
#     }
#     write.table(head,file=paste(root,dest,"/pred.csv",sep=""),col.names=FALSE,quote=FALSE,row.names=FALSE)
#     IDS <- rownames(p)
#     p <- cbind(IDS,p)
#     write.table(p,file=paste(root,dest,"/pred.csv",sep=""),quote=FALSE,row.names=F,sep=",",append=TRUE)
#
#     #Create combined probabilities
#     if (classify=="Y") {
#       head<-","
#       p<-data.frame()
#       for (i in 1:length(titles)) {
#         head<-paste(head,",",titles[i],sep="")
#         if (ncol(probs[[i]])>2)
#           for (j in 2:(ncol(probs[[i]])-1))
#             head<-paste(head,",",sep="")
#           if (ncol(p)==0) p<-data.frame(probs[[i]])
#           else p<-data.frame(p,subset(probs[[i]],select=-1),check.names=FALSE)
#       }
#       write.table(head,file=paste(root,dest,"/prob.csv",sep=""),col.names=FALSE,quote=FALSE,row.names=FALSE)
#       IDS <- rownames(p)
#       p <- cbind(IDS,p)
#       write.table(p,file=paste(root,dest,"/prob.csv",sep=""),quote=FALSE,row.names=FALSE,sep=",",append=TRUE)
#     }


#     #Create combined summaries
#     sums<-c()
#     if (bndir!="") sums<-c(sums,readLines(paste(root,bndir,"/summary.txt",sep="")))
#     if (phdir!="") sums<-c(sums,readLines(paste(root,phdir,"/summary.txt",sep="")))
#     if (apdir!="") sums<-c(sums,readLines(paste(root,apdir,"/summary.txt",sep="")))
#     if (fpdir!="") sums<-c(sums,readLines(paste(root,fpdir,"/summary.txt",sep="")))
#     if (chdir!="") sums<-c(sums,readLines(paste(root,chdir,"/summary.txt",sep="")))
#     writeLines(sums,paste(root,dest,"/summary.txt",sep=""))
  # }
}


SplitAnova <- function(splitdata, metric, file = NA) {

  single.desc <- (length(unique(splitdata$Descriptor)) == 1)
  # Calculate anova
  out <- glm(Model.Acc ~ Split + Trmt, data = splitdata)
  # total df - df for residuals
  mod.df <- out$df.null - out$df.residual
  # is this TSS - SSR?
  mod.dev <- out$null.deviance - out$deviance
  mod.ms <- mod.dev/mod.df
  err.df <- out$df.residual
  err.dev <- out$deviance
  err.ms <- err.dev/err.df
  tot.df <- out$df.null
  tot.dev <- out$null.deviance
  f.stat <- mod.ms/err.ms
  p.val <- df(f.stat, mod.df, err.df)
  if (p.val < 1e-04)
    p.val <- "<.0001"
  form.df <- format(c("DF", mod.df, err.df, tot.df), width = 3, justify = "right")
  form.dec.num <- format(c(mod.dev, err.dev, tot.dev, mod.ms, err.ms, f.stat),
    digits = 4, nsmall = 3)
  form.dec <- format(c("SS", "MS", "F", form.dec.num), digits = 4, nsmall = 3,
    justify = "right")
  form.p <- format(c("p-value", if (is.numeric(p.val)) round(p.val, digits = 4) else p.val),
    digits = 1, nsmall = 4, justify = "right")
  cat(paste0("   Analysis of Variance on: '",metric,"'\n"))
  if (single.desc)
    cat(paste("            Using factors: Split and Method\n"))
  else cat(paste(" Using factors: Split and Descriptor/Method combination\n"))
  cat(paste("Source", form.df[1], form.dec[1], form.dec[2], form.dec[3], form.p[1],
    "\n", sep = "   "))
  cat(paste("Model ", form.df[2], form.dec[4], form.dec[7], form.dec[9], form.p[2],
    "\n", sep = "   "))
  cat(paste("Error ", form.df[3], form.dec[5], form.dec[8], "\n", sep = "   "))
  cat(paste("Total ", form.df[4], form.dec[6], "\n", sep = "   "))

  if (!is.na(file))
    pdf(file)

  # plot(c(-1,1),c(-1,1),xlab='',ylab='',axes=FALSE,type='n')
  # text(-1,.9,str1,adj=c(0,.5),family='mono',cex=.7)
  # text(-1,.8,str2,adj=c(0,.5),family='mono',cex=.7)
  # text(rep(-1,4),seq(.6,.3,by=-.1),c(str3,str4,str5,str6),adj=c(0,0),family='mono',cex=.7)
  root.mse <- sqrt(err.ms)
  all.mean <- mean(splitdata$Model.Acc)
  coef.var <- root.mse/all.mean * 100
  r.sq <- mod.dev/tot.dev
  form.stat.num <- format(c(r.sq, coef.var, root.mse, all.mean), digits = 3, nsmall = 4)
  form.stat <- format(c("R-Square", "Coef Var", "Root MSE", "Mean", form.stat.num),
    digits = 3, nsmall = 4, justify = "right")
  cat(paste("   ", form.stat[1], form.stat[2], form.stat[3], form.stat[4], "\n",
    sep = "   "))
  cat(paste("   ", form.stat[5], form.stat[6], form.stat[7], form.stat[8], "\n",
    sep = "   "))
  # text(rep(-1,2),seq(.1,0,by=-.1),c(str7,str8),adj=c(0,0),family='mono',cex=.7)

  aout <- anova(out)
  f.df <- aout$Df[2]
  f.dev <- aout$Deviance[2]
  f.ms <- f.dev/f.df
  f.f.stat <- f.ms/err.ms
  t.df <- aout$Df[3]
  t.dev <- aout$Deviance[3]
  t.ms <- t.dev/t.df
  t.f.stat <- t.ms/err.ms
  f.p.val <- df(f.f.stat, f.df, t.df)
  if (f.p.val < 1e-04)
    f.p.val <- "<.0001"
  t.p.val <- df(t.f.stat, t.df, err.df)
  if (t.p.val < 1e-04)
    t.p.val <- "<.0001"
  form.df <- format(c("DF", f.df, t.df), width = 3, justify = "right")
  form.dec.num <- format(c(f.dev, f.ms, f.f.stat, t.dev, t.ms, t.f.stat), digits = 3,
    nsmall = 3)
  form.dec <- format(c("SS", "MS", "F", form.dec.num), digits = 3, nsmall = 3,
    justify = "right")
  form.p <- format(c("p-value", f.p.val, t.p.val), digits = 1, nsmall = 4, justify = "right")
  form.p <- format(c("p-value", if (is.numeric(f.p.val)) round(f.p.val, digits = 4) else f.p.val,
    if (is.numeric(t.p.val)) round(t.p.val, digits = 4)
    else t.p.val), digits = 1, nsmall = 4, justify = "right")
  cat(paste("Source   ", form.df[1], form.dec[1], form.dec[2], form.dec[3], form.p[1],
    "\n", sep = "   "))
  cat(paste("Split    ", form.df[2], form.dec[4], form.dec[5], form.dec[6], form.p[2],
    "\n", sep = "   "))
  if (single.desc)
    cat(paste("Method   ", form.df[3], form.dec[7], form.dec[8], form.dec[9],
      form.p[3], "\n", sep = "   "))
  else cat(paste("Desc/Meth", form.df[3], form.dec[7], form.dec[8], form.dec[9],
    form.p[3], "\n", sep = "   "))
  # text(rep(-1,3),seq(-.2,-.4,by=-.1),c(str9,str10,str11),adj=c(0,0),family='mono',cex=.7)

  # Calculate lsmeans
  lsmeans <- c()
  fit <- fitted.values(out)
  for (i in levels(splitdata$Trmt)) lsmeans <- c(lsmeans, mean(fit[splitdata$Trmt ==
    i]))

  # Calculate Tukey-Kramer adjusted p-values for diff means
  aov.out <- aov(Model.Acc ~ Split + Trmt, data = splitdata)
  tout <- TukeyHSD(aov.out, "Trmt")

  # Create matrix of p-values
  n <- length(levels(splitdata$Trmt))
  pval <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      pval[i, j] <- tout$Trmt[n * (i - 1) - i * (i - 1)/2 + j - i, 4]
      pval[j, i] <- pval[i, j]
    }
  }

  # Plot
  McsPlot(lsmeans, pval, as.numeric(levels(splitdata$Trmt)),
    as.character(unique(splitdata$Descriptor)),
    as.character(unique(splitdata$Method)), single.desc, metric)
  if (!is.na(file))
    dev.off()
}

CombineSplits <- function(cml.result, file = NA, metric = "enhancement", at = NA, thresh = 0.5, ...) {
  y <- cml.result$responses
  # take initial enhancement at 300 tests or if less then 300, the number of y's/4
  if (is.na(at)) {
    at <- min(300, ceiling(length(y)/4))
  } else if (at > length(y)) {
    stop("'at' needs to be smaller than the number of responses")
  }

  # makes desciptor set names shorter so that they fit on the MCS plot
  abbrev.names <- c()
  num.desc <- length(cml.result$all.preds[[1]])
  des.names <- names(cml.result$all.preds[[1]])
  if (grepl("Descriptor Set", des.names[1])) {
    for (i in 1:num.desc) {
      # TO DO add option to specify abbreviated names?
      abbrev.names <- c(abbrev.names, paste0("Des", i))
    }
  } else {
    for (i in 1:num.desc) {
      abbrev.names <- c(abbrev.names, substr(des.names[i], 1, 4))
    }
  }
  if (cml.result$classify == "N") {
    # then no classification model
    # TODO: Why are they concerned about these objects
    # existing?
    if (exists("sp"))
      rm(sp)
    if (exists("ds"))
      rm(ds)
    if (exists("me"))
      rm(me)
    if (exists("ma"))
      rm(ma)

    # TODO Check model accuracy measures
    for (split in 1:length(cml.result$all.preds)) {
      pred <- cml.result$all.preds[[split]]
      for (i in 1:length(pred)) {
        desc <- abbrev.names[i]
        for (j in 2:ncol(pred[[i]])) {
          if (metric == "enhancement") {
            model.acc <- EnhancementCont(pred[[i]][, j], y, at)
          } else if (metric == "R2") {
            model.acc <- (cor(y, pred[[i]][, j]))^2
          } else if (metric == "RMSE") {
            model.acc <- sqrt(sum(y - pred[[i]][, j])^2)
          } else if (metric == "rho") {
            model.acc <- cor(y, pred[[i]][, j], method = "spearman")
          } else {
            stop("y is continuous. 'metric' should be a model accuracy measure
              implemented for continuous response in ChemModLab")
          }
          colnm <- names(pred[[i]])[j]
          # TO DO do we still need this?
          period <- regexpr(".", colnm, fixed = TRUE)[1]
          if (period > 0)
            meth <- substr(colnm, 1, (period - 1))
          else meth <- colnm
          if (!exists("sp"))
            sp <- split
          else sp <- c(sp, split)
          if (!exists("ds"))
            ds <- desc
          else ds <- c(ds, desc)
          if (!exists("me"))
            me <- meth
          else me <- c(me, meth)
          if (!exists("ma"))
            ma <- model.acc
          else ma <- c(ma, model.acc)
          }
      }
    }
  } else {
    if (exists("sp"))
      rm(sp)
    if (exists("ds"))
      rm(ds)
    if (exists("me"))
      rm(me)
    if (exists("ma"))
      rm(ma)
    for (split in 1:length(cml.result$all.preds)) {
      prob <- cml.result$all.probs[[split]]
      pred <- cml.result$all.preds[[split]]
      for (i in 1:length(prob)) {
        desc <- abbrev.names[i]
        for (j in 2:ncol(prob[[i]])) {
          # TODO: Determine whether to use predicted probabilities or predictions
          # for models that provide both - sometimes they differ.
          yhat <- prob[[i]][, j] > thresh
          if (metric == "enhancement") {
            model.acc <- Enhancement(prob[[i]][, j], y, at)
          } else if (metric == "auc") {
            model.acc <- as.numeric(auc(y, prob[[i]][, j]))
          } else if (metric == "error rate") {
            model.acc <- mean(y != yhat)
          } else if (metric == "specificity") {
            idx <- y == 0
            model.acc <- mean(y[idx] == yhat[idx])
          } else if (metric == "sensitivity") {
            idx <- y == 1
            model.acc <- mean(y[idx] == yhat[idx])
            # TODO: Remove rho?
#           } else if (metric == "rho") {
#             model.acc <- cor(y, prob[[i]][, j], method = "spearman")
          } else {
            stop("y is binary. 'metric' should be a model accuracy measure
              implemented for binary response in ChemModLab")
          }
          colnm <- names(prob[[i]])[j]
          period <- regexpr(".", colnm, fixed = TRUE)[1]
          if (period > 0)
            meth <- substr(colnm, 1, (period - 1))
          else meth <- colnm
          if (!exists("sp"))
            sp <- split else sp <- c(sp, split)
          if (!exists("ds"))
            ds <- desc else ds <- c(ds, desc)
          if (!exists("me"))
            me <- meth else me <- c(me, meth)
          if (!exists("ma"))
            ma <- model.acc else ma <- c(ma, model.acc)
          }
      }
      for (i in 1:length(pred)) {
        desc <- abbrev.names[i]
        for (j in 2:ncol(pred[[i]])) {
          yhat <- pred[[i]][, j] > thresh
          if (metric == "enhancement") {
            model.acc <- Enhancement(pred[[i]][, j], y, at)
          } else if (metric == "auc") {
            model.acc <- as.numeric(auc(y, pred[[i]][, j]))
          } else if (metric == "error rate") {
            model.acc <- mean(y != yhat)
          } else if (metric == "specificity") {
            idx <- y == 0
            model.acc <- mean(y[idx] == yhat[idx])
          } else if (metric == "sensitivity") {
            idx <- y == 1
            model.acc <- mean(y[idx] == yhat[idx])
          } else if (metric == "rho") {
            model.acc <- cor(y, pred[[i]][, j], method = "spearman")
          } else {
            stop("y is binary. 'metric' should be a model accuracy measure
              implemented for binary response in ChemModLab")
          }
          colnm <- names(pred[[i]])[j]
          # TODO: do we still need this?
          period <- regexpr(".", colnm, fixed = TRUE)[1]
          if (period > 0)
            meth <- substr(colnm, 1, (period - 1)) else meth <- colnm
          if (!exists("sp"))
            sp <- split else sp <- c(sp, split)
          if (!exists("ds"))
            ds <- desc else ds <- c(ds, desc)
          if (!exists("me"))
            me <- meth else me <- c(me, meth)
          if (!exists("ma"))
            ma <- model.acc else ma <- c(ma, model.acc)
          }
        }
      }
    }

  out <- data.frame(sp, ds, me, ma)
  names(out) <- c("Split", "Descriptor", "Method", "Model.Acc")

  # for performance reasons it would be preferable not to calculate the enhancement
  # for the classification methods, but for speed of implementation this is being
  # done here.  Dropping all the model accuracy calculated using predictions from
  # classification methods

  out <- out[!duplicated(out[, 1:3]), ]
  out$Split <- factor(out$Split)

  methods <- as.character(unique(out$Method))
  num.enhs <- dim(out)[1]
  descriptors <- abbrev.names
  out$Trmt <- vector(mode = "logical", length = num.enhs)
  for (i in 1:num.enhs) {
    out$Trmt[i] <- which(methods == out$Method[i])
  }
  for (i in 1:num.enhs) {
    out$Trmt[i] <- out$Trmt[i] + 100 * (which(descriptors == out$Descriptor[i]))
  }

  #   if (length(unique(out$Descriptor))==1){
  #     out$Trmt<-1*(out$Method=="ENet") + 2*(out$Method=="RF") + 3*(out$Method=="KNN") + 4*(out$Method=="LAR") +
  #     5*(out$Method=="NNet") + 6*(out$Method=="PCR") + 7*(out$Method=="PLS") + 8*(out$Method=="PLSLDA") +
  #     9*(out$Method=="RPart") + 10*(out$Method=="Ridge") + 11*(out$Method=="SVM") + 12*(out$Method=="Tree")
  #   } else{
  #     out$Trmt<-100*( 1*(out$Descriptor=="Atom Pairs") + 2*(out$Descriptor=="Burden Numbers") + 3*(out$Descriptor=="Carhart Atom Pairs")
  #                     + 4*(out$Descriptor=="Fragment Pairs") + 5*(out$Descriptor=="Pharmacophores") ) +
  #     1*(out$Method=="ENet") + 2*(out$Method=="RF") + 3*(out$Method=="KNN") + 4*(out$Method=="LAR") +
  #     5*(out$Method=="NNet") + 6*(out$Method=="PCR") + 7*(out$Method=="PLS") + 8*(out$Method=="PLSLDA") +
  #     9*(out$Method=="RPart") + 10*(out$Method=="Ridge") + 11*(out$Method=="SVM") + 12*(out$Method=="Tree")
  #   }

  out$Trmt <- factor(out$Trmt)
  SplitAnova(out, metric, file)
}

# Suspending heatmaps at the moment-------------------

#
# # --- Redefining the R heatmap function, to handle our tweaks
# heatmap<-function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL,
#                    distfun = dist, hclustfun = hclust, reorderfun = function(d,
#                                                                              w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv,
#                                                                                                                                         "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE,
#                    margins = c(1, 5), ColSideColors, RowSideColors, cexRow = 1/log10(nr),
#                    cexCol = 1, labRow = NULL,
#                    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE,
#                    verbose = getOption("verbose"), ...)
# {
#   scale <- if (symm && missing(scale))
#     "none"
#   else match.arg(scale)
#   if (length(di <- dim(x)) != 2 || !is.numeric(x))
#     stop("'x' must be a numeric matrix")
#   nr <- di[1]
#   nc <- di[2]
#   if (nr <= 1 || nc <= 1)
#     stop("'x' must have at least 2 rows and 2 columns")
#   if (!is.numeric(margins) || length(margins) != 2)
#     stop("'margins' must be a numeric vector of length 2")
#   doRdend <- !identical(Rowv, NA)
#   doCdend <- !identical(Colv, NA)
#   if (is.null(Rowv))
#     Rowv <- rowMeans(x, na.rm = na.rm)
#   if (is.null(Colv))
#     Colv <- colMeans(x, na.rm = na.rm)
#   if (doRdend) {
#     if (inherits(Rowv, "dendrogram"))
#       ddr <- Rowv
#     else {
#       hcr <- hclustfun(distfun(x))
#       ddr <- as.dendrogram(hcr)
#       if (!is.logical(Rowv) || Rowv)
#         ddr <- reorderfun(ddr, Rowv)
#     }
#     if (nr != length(rowInd <- order.dendrogram(ddr)))
#       stop("row dendrogram ordering gave index of wrong length")
#   }
#   else rowInd <- 1:nr
#   if (doCdend) {
#     if (inherits(Colv, "dendrogram"))
#       ddc <- Colv
#     else if (identical(Colv, "Rowv")) {
#       if (nr != nc)
#         stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
#       ddc <- ddr
#     }
#     else {
#       hcc <- hclustfun(distfun(if (symm)
#         x
#         else t(x)))
#       ddc <- as.dendrogram(hcc)
#       if (!is.logical(Colv) || Colv)
#         ddc <- reorderfun(ddc, Colv)
#     }
#     if (nc != length(colInd <- order.dendrogram(ddc)))
#       stop("column dendrogram ordering gave index of wrong length")
#   }
#   else colInd <- 1:nc
#   x <- x[rowInd, colInd]
#   labRow <- if (is.null(labRow))
#     if (is.null(rownames(x)))
#       (1:nr)[rowInd]
#   else rownames(x)
#   else labRow[rowInd]
#   labCol <- if (is.null(labCol))
#     if (is.null(colnames(x)))
#       (1:nc)[colInd]
#   else colnames(x)
#   else labCol[colInd]
#   if (scale == "row") {
#     x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
#     sx <- apply(x, 1, sd, na.rm = na.rm)
#     x <- sweep(x, 1, sx, "/")
#   }
#   else if (scale == "column") {
#     x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
#     sx <- apply(x, 2, sd, na.rm = na.rm)
#     x <- sweep(x, 2, sx, "/")
#   }
#   lmat <- rbind(c(NA, 3), 2:1)
#   lwid <- c(if (doRdend) 1 else 0.05, 4)
#   lhei <- c(if (doCdend) 1 else 0.25, 4)
#   if (!missing(ColSideColors)) {
#     if (!is.character(ColSideColors) || length(ColSideColors) !=
#         nc)
#       stop("'ColSideColors' must be a character vector of length ncol(x)")
#     lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
#     lhei <- c(lhei[1], 0.2, lhei[2])
#   }
#   if (!missing(RowSideColors)) {
#     if (!is.character(RowSideColors) || length(RowSideColors) !=
#         nr)
#       stop("'RowSideColors' must be a character vector of length nrow(x)")
#     lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1),
#                                    1), lmat[, 2] + 1)
#     lwid <- c(lwid[1], 0.2, lwid[2])
#   }
#   lmat[is.na(lmat)] <- 0
#   if (verbose) {
#     cat("layout: widths = ", lwid, ", heights = ", lhei,
#         "; lmat=\n")
#     print(lmat)
#   }
#   op <- par(no.readonly = TRUE)
#   on.exit(par(op))
#   layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
#   if (!missing(RowSideColors)) {
#     par(mar = c(margins[1], 0, 0, 0.5))
#     image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
#   }
#   if (!missing(ColSideColors)) {
#     par(mar = c(0.5, 0, 0, margins[2]))
#     image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
#   }
#   par(mar = c(margins[1], 0, 0, margins[2]))
#   if (!symm || scale != "none")
#     x <- t(x)
#   if (revC) {
#     iy <- nr:1
#     ddr <- rev(ddr)
#     x <- x[, iy]
#   }
#   else iy <- 1:nr
#   image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
#           c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
#   axis(1, 1:nc, labels = labCol, las = 2, line = -0.4, tick = 0,
#        cex.axis = cexCol)
#   if (!is.null(xlab))
#     mtext(xlab, side = 1, line = margins[1] - 1.25)
#   axis(4, iy[seq(1,length(labRow),by=2)], labels = labRow[seq(1,length(labRow),by=2)],
#        las = 2, line = -0.4, tick = 0, cex.axis = cexRow)
#   axis(4, iy[seq(2,length(labRow),by=2)], labels = labRow[seq(2,length(labRow),by=2)],
#        las = 2, line = 1, tick = 0, cex.axis = cexRow)
#   if (!is.null(ylab))
#     mtext(ylab, side = 4, line = margins[2] - 1.25)
#   if (!missing(add.expr))
#     eval(substitute(add.expr))
#   par(mar = c(margins[1], 0, 0, 0))
#   if (doRdend)
#     plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
#   else frame()
#   par(mar = c(0, 0, 1, margins[2]))
#   if (doCdend)
#     plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
#   else if (!is.null(main))
#     frame()
#   if (!is.null(main))
#     title(main, cex.main = 1.5)
#   invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro &&
#                                                               doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
# }


# Disregarding prediction at the moment

BackPredict <- function(yfilein, ycol, xfilein, xcols, filepred, fileimpdesc, fileprob,
  infofile = NA, idcol = NA, modelyfile, modelxfile) {
  #-----Background predicting program
  # yfilein=input response data: response in column ycol xfilein=input descriptors
  # data: descriptors in columns xcols filepred=output file of predictions:
  # observed response in column 1, one additional column of predictions for each
  # method modelyfile=model response data modelxfile=model descriptors data
  # idcol=column which contains compound labels

  #-----Read in methods
  if (missing(infofile))
    meths <- c("NNet", "PCR", "ENet", "PLS", "Ridge", "LARs", "PLSLDA", "RPart",
      "Tree", "SVM", "KNN", "Forest", "Forest70", "TreeEns", "RPartEns", "KNNEns")
  else meths <- read.table(file = infofile, sep = "|", skip = 1, nrows = 1, colClasses = "character")

  #-----Read in model data
  model.list <- read.in.data(yfilein = modelyfile, xfilein = modelxfile, idcol = idcol)
  model.data <- model.list[[1]]

  #-----Read in prediction data, and discard constant columns from model data
  work.data <- read.in.data(yfilein, ycol, xfilein, xcols, idcol = idcol, type = "PREDICT")[[1]]
  work.data <- subset(work.data, select = (!model.list[[2]]))

  n.obs <- nrow(work.data)
  n.pred <- ncol(work.data) - 1



  #-----Determine if we have binary or continuous data
  if (sum(!(model.data$y %in% c(1, 0))))
    classify <- "N" else classify <- "Y"

  #-----Start outputing to summary file
  cat("Begining Prdeictions: ", date(), "\n")
  cat("Model Descriptors: ", modelxfile, "(", (ncol(model.data) - 1), ")\n")
  cat("Model Responses: ", modelyfile, "(", nrow(model.data), ")\n")
  cat("Prediction Descriptors: ", xfilein, "(", n.pred, ")\n")
  if ((!missing(yfilein)) && (yfilein != ""))
    cat("Prediction Responses: ", yfilein, "(", n.obs, ")\n\n\n")
  all.preds <- data.frame(Observed = work.data$y)
  if (classify == "Y")
    all.probs <- data.frame(Observed = work.data$y) else all.probs <- NA

  IDS <- rownames(work.data)

  #-----'tree' method
  if (sum(meths %in% "Tree") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredTree(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...Tree not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, Tree = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, Tree = work.results$prob)
    }
  }

  #-----'rpart' method
  if (sum(meths %in% "RPart") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredRpartCml(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...RPart not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, RPart = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, RPart = work.results$prob)
    }
  }

  #-----'randomforest' method
  if (sum(meths %in% "Forest") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredRf(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...RF not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, RF = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, RF = work.results$prob)
    }
  }

  #-----'svm' method
  if (sum(meths %in% "SVM") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredSvm(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...SVM not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, SVM = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, SVM = work.results$prob)
    }
  }

  #-----'nnet' method
  if ((sum(meths %in% "NNet") == 1) && (classify == "Y")) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredNnet(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...NNet not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, NNet = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, NNet = work.results$prob)
    }
  }

  #-----'knn' method
  if ((sum(meths %in% "KNN") == 1) && (classify == "Y")) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredKnn(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...KNN not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, KNN = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, KNN = work.results$prob)
    }
  }

  #-----'pls.lda' method
  if ((sum(meths %in% "PLSLDA") == 1) && (classify == "Y")) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredPlsLdaNew(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...PLSLDA not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, PLSLDA = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, PLSLDA = work.results$prob)
    }
  }

  #-----'lars' method
  if (sum(meths %in% "LARs") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredLars(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...LAR not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, LAR = work.results$pred)
    }
  }

  #-----'lm.ridge' method
  if (sum(meths %in% "Ridge") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredRidge(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...Ridge not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, Ridge = work.results$pred)
    }
  }

  #-----'enet' method
  if (sum(meths %in% "ENet") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredEnet(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...ENet not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, ENet = work.results$pred)
    }
  }

  #-----'PcrZG' method
  if (sum(meths %in% "PCR") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredPcr(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...PCR not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, PCR = work.results$pred)
    }
  }

  #-----'pls.R' method
  if (sum(meths %in% "PLS") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredPlsR(model.data, work.data,
      classify), error = function(e) {
      warning(paste("WARNING...PLS not run:", e$message))
      work.results <- list()
    }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, PLS = work.results$pred)
    }
  }

  cat("Ending Predictions: ", date(), "\n\n\n")
  #-----Write 'predictions', 'probabilities' and 'important descriptors' files
  desc <- xfilein
  while ((pos <- regexpr("/", desc, fixed = TRUE)) > 0) desc <- substr(desc, (pos +
    1), nchar(desc))
  desc <- substr(desc, 1, (regexpr(".", desc, fixed = TRUE) - 1))
  head <- paste(",,", desc, sep = "")
  if (ncol(all.preds) > 2)
    for (i in 2:(ncol(all.preds) - 1)) head <- paste(head, ",", sep = "")
  write.table(head, file = filepred, quote = FALSE, col.names = FALSE, row.names = FALSE)
  all.preds <- cbind(IDS, all.preds)
  write.table(all.preds, file = filepred, quote = FALSE, sep = ",", row.names = FALSE,
    append = TRUE)
  if (classify == "Y") {
    head <- paste(",,", desc, sep = "")
    if (ncol(all.probs) > 2)
      for (i in 2:(ncol(all.probs) - 1)) head <- paste(head, ",", sep = "")
    write.table(head, file = fileprob, quote = FALSE, col.names = FALSE, row.names = FALSE)
    all.probs <- cbind(IDS, all.probs)
    write.table(all.probs, file = fileprob, quote = FALSE, sep = ",", row.names = FALSE,
      append = TRUE)
  }

  return(list(all.preds = all.preds, all.probs = all.probs, classify = classify,
    responses = work.data$y))
}


PredTree <- function(model.data, work.data, classify) {
  #-----Recursive partitioning using 'tree' with splitting criterion deviance and default settings. Specifically:
  # mincut=5, minimum leaf size minsize=10, minimum parent size mindev=.01,
  # within-node deviance must be at least this times that of the root node for node
  # to split
  cat("Tree ---------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- tree(as.factor(y) ~ ., data = model.data, method = "class")
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
      work.data, type = "class")]
    work.prob <- predict(work.meth, work.data)[, 1]
  } else {
    work.meth <- tree(y ~ ., data = work.data)
    work.pred <- predict(work.meth, work.data, type = "vector")
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}

#--- RPart has an internal function named 'PredRpart', our function is PredRpart.chm
PredRpartCml <- function(model.data, work.data, classify) {
  #-----Recursive partitioning using 'rpart' with splitting criterion 'information' and
  # minbucket=5, minimum leaf size minsplit=10, minimum parent size maxcompete=0,
  # don't get information on competitive splits maxsurrogate=0, don't get
  # information on surrogate splits Possible modifications that have NOT been
  # pursued here: many ...
  cat("RPart --------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- rpart(as.factor(y) ~ ., data = model.data, method = "class",
      parms = list(split = "information"), minsplit = 10, minbucket = 5, maxcompete = 0,
      maxsurrogate = 0)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
      work.data, type = "class")]
    work.prob <- predict(work.meth, work.data)[, 1]
  } else {
    work.meth <- rpart(y ~ ., data = model.data, method = "anova", parms = list(split = "information"),
      minsplit = 10, minbucket = 5, maxcompete = 0, maxsurrogate = 0)
    work.pred <- predict(work.meth, work.data, type = "vector")
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredRf <- function(model.data, work.data, classify) {
  #-----Random Forest using
  # ntree=100 mtry = sqrt(p) [for classification] and = p/3 [for regression]
  # nodesize = 5 importance = TRUE (to calculate important descriptors) Possible
  # modifications that have NOT been pursued here: many ...
  cat("Forest -------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- randomForest(y = as.factor(model.data$y), x = model.data[, -1],
      ntree = 100, nodesize = 5)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
      work.data[, -1])]
    work.prob <- predict(work.meth, work.data[, -1], type = "prob")[, 1]
  } else {
    work.meth <- randomForest(y = model.data$y, x = model.data[, -1], ntree = 100,
      nodesize = 5)
    work.pred <- predict(work.meth, work.data[, -1])
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredSvm <- function(model.data, work.data, classify) {
  #-----Support Vector Machine using
  # kernel = radial.basis gamma = 1 Possible modifications that have NOT been
  # pursued here: many ...
  cat("SVM ----------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- svm(y = as.factor(model.data$y), x = model.data[, -1], gamma = 1,
      probability = TRUE)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
      work.data[, -1])]
    temp.prob <- predict(work.meth, work.data[, -1], probability = TRUE)
    work.prob <- attr(temp.prob, "probabilities")[, 1]
  } else {
    work.meth <- svm(y = model.data$y, x = model.data[, -1], gamma = 1)
    work.pred <- predict(work.meth, work.data[, -1])
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredNnet <- function(model.data, work.data, classify) {
  #-----Neural Network using 'nnet' with:
  # size=2 trace=FALSE (don't print out convergence info) Possible modifications
  # that have NOT been pursued here: many
  cat("NNet ---------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.pred <- rep(NA, n.work)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- nnet(as.factor(y) ~ ., data = model.data, size = 2, trace = FALSE)
    work.pred <- predict(work.meth, work.data, type = "class")
    work.prob <- 1 - predict(work.meth, work.data, type = "raw")
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredKnn <- function(model.data, work.data, classify) {
  #-----K nearest neighbor using
  # k = 3 Possible modifications that have NOT been pursued here: cv option ...
  cat("KNN ----------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.pred <- rep(NA, n.work)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- knn(cl = as.factor(model.data$y), train = model.data[, -1],
      test = work.data[, -1], k = 10, prob = TRUE)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[work.meth]
    work.prob <- abs(work.pred + attr(work.meth, "prob") - 1)
  }
  return(list(pred = work.pred, prob = work.prob))
}


PredKnnflex <- function(model.data, work.data, classify) {
  #-----K nearest neighbor using
  # k = 3
  cat("KNN ----------------------------------------", "\n")
  work.meth <- knn.dist(rbind(model.data[, -1], work.data[, -1]))
  n.model <- nrow(model.data)
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.pred <- as.numeric(knn.predict(train = 1:n.model, test = (n.model +
      1):(n.work + n.model), y = work.data$y, k = 10, dist.matrix = work.meth,
      agg.meth = "majority"))
    work.prob <- knn.predict(train = 1:n.model, test = (n.model + 1):(n.work +
      n.model), y = work.data$y, k = 10, dist.matrix = work.meth, agg.meth = "mean")
  } else {
    work.pred <- knn.predict(train = 1:n.model, test = (n.model + 1):(n.work +
      n.model), y = work.data$y, k = 10, dist.matrix = work.meth, agg.meth = "mean")
  }
  return(list(pred = work.pred, prob = work.prob))
}


PredLars <- function(model.data, work.data, classify) {
  #-----Lars
  cat("LARs ---------------------------------------", "\n")
  work.meth <- lars(y = model.data[, 1], x = as.matrix(model.data[, -1]), type = "lar")
  temp.pred <- predict(work.meth, as.matrix(work.data[, -1]), type = "fit")
  work.pred <- temp.pred$fit[, dim(temp.pred$fit)[2]]
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredRidge <- function(model.data, work.data, classify) {
  #-----Ridge Regression using
  # lambda = 0.1
  cat("Ridge --------------------------------------", "\n")
  work.meth <- lm.ridge(y ~ ., data = model.data, lambda = 0.1)
  work.pred <- (as.matrix(work.data[, -1]) %*% (work.meth$coef/work.meth$scales)) +
    ((work.meth$ym - (work.meth$xm %*% (work.meth$coef/work.meth$scales)))[1,
      1])
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredEnet <- function(model.data, work.data, classify) {
  #-----enet
  cat("ENet ---------------------------------------", "\n")
  work.meth <- enet(y = model.data$y, x = as.matrix(model.data[, -1]), lambda = 1)
  temp.pred <- predict(work.meth, as.matrix(work.data[, -1]), type = "fit")
  work.pred <- temp.pred$fit[, dim(temp.pred$fit)[2]]
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredPcr <- function(model.data, work.data, classify) {
  #-----Principal components regression using home-grown code
  cat("PCR ----------------------------------------", "\n")
  newX <- model.data[, -1]
  work.meth <- PcrZG(X = model.data[, -1], Y = model.data$y, newX = work.data[,
    -1])
  work.pred <- work.meth$Ypred
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredPlsR <- function(model.data, work.data, classify) {
  #-----Partial least squares using home-grown code based on 'kernelpls'
  cat("PLS ----------------------------------------", "\n")
  work.meth <- KernelPlsNew(X = model.data[, -1], Y = model.data$y, ncomp = min(nrow(model.data),
    (ncol(model.data) - 1), 100), newX = work.data[, -1])
  nLV.ZG <- ZhuGhodsi(work.meth$Yvar)
  work.pred <- work.meth$Ypred[, 1, nLV.ZG]
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredPlsLdaNew <- function(model.data, work.data, classify) {
  #-----Partial least squares using home-grown code based on 'kernelpls'
  cat("PLSLDA -------------------------------------", "\n")
  work.meth <- SimPlsLdaNew(X = model.data[, -1], Y = model.data$y, ncomp = min(nrow(model.data),
    (ncol(model.data) - 1), 100), newX = work.data[, -1])
  nLV.ZG <- max(ZhuGhodsi(work.meth$Yvar), 1)
  work.prob <- work.meth$Yprob[, 2, nLV.ZG]
  work.pred <- work.meth$Ypred[, nLV.ZG]
  return(list(pred = work.pred, prob = work.prob))
}


