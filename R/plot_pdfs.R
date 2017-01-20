plotBB <- function(bb, xfilein, plotfile) {
  y <- bb$responses[, 1]
  desc <- read.csv(xfilein, row.names = 1)
  # bitmap(paste(plotfile,'.bmp',sep=''))
  pdf(paste(plotfile, ".pdf", sep = ""))
  if ((bb$classify == "Y") && (ncol(bb$all.probs) > 0)) {
    HitCurve(bb$all.probs[, -1], y = y)
    ContCurve(bb$all.preds[, !(names(bb$all.preds) %in% names(bb$all.probs))],
              y = y, curves.only = TRUE, start.col = (ncol(bb$all.probs) - 1))
  } else {
    ContCurve(bb$all.preds[, -1], y = y)
  }
  dev.off()
  pdf(paste(plotfile, "s.pdf", sep = ""))
  if ((bb$classify == "Y") && (ncol(bb$all.probs) > 0)) {
    HitCurve(bb$all.probs[, -1], y = y)
    ContCurve(bb$all.preds[, !(names(bb$all.preds) %in% names(bb$all.probs))],
              y = y, curves.only = TRUE, start.col = (ncol(bb$all.probs) - 1))
    HitMap(data.frame(bb$all.probs[, -1], bb$all.preds[, !(names(bb$all.preds) %in%
                                                             names(bb$all.probs))]), y = y, y.labels = row.names(desc), descriptors = desc)
  } else {
    ContCurve(bb$all.preds[, -1], y = y)
    HitMap(bb$all.preds[, -1], y = y, y.labels = row.names(desc), descriptors = desc)
  }
  dev.off()
}
