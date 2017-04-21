HitRate <- function(prob, y, at) {
  # hit rate is accumulation/n for each pt
  hrpts <- Accumulation(prob, y)/(1:length(y))
  
  # return the hit rate for every point, or those specified by 'at'
  if (missing(at))
    return(hrpts) else return(hrpts[at])
}

Enhancement <- function(prob, y, at) {
  # enhancement is the hit rate/(M/N)
  epts <- HitRate(prob, y)/(sum(y)/length(y))
  
  # return the enhancement for every point, or those specified by 'at'
  if (missing(at))
    return(epts) else return(epts[at])
}

EnhancementCont <- function(pred, y, at) {
  # assumes bigger y is better to avoid problems with negative numbers, convert y
  # to y-min(y) 'at' assumed to be single number, not vector enhancement@at is
  # (mean y over top at)/(mean y over all)
  y <- y - min(y)
  pred.order <- order(pred, decreasing = TRUE)
  enh <- mean(y[pred.order[1:at]])/mean(y)
  return(enh)
}

BackSpecificity <- function(pred, yhat, y, at) {
  pred.order <- order(pred, decreasing = TRUE)
  y <- y[pred.order[1:at]]
  yhat <- yhat[pred.order[1:at]]
  idx <- y == 0
  mean(y[idx] == yhat[idx])
}

BackSensitivity <- function(pred, yhat, y, at) {
  pred.order <- order(pred, decreasing = TRUE)
  y <- y[pred.order[1:at]]
  yhat <- yhat[pred.order[1:at]]
  idx <- y == 1
  mean(y[idx] == yhat[idx])
}

BackSpecificity <- function(pred, yhat, y, at) {
  pred.order <- order(pred, decreasing = TRUE)
  y <- y[pred.order[1:at]]
  yhat <- yhat[pred.order[1:at]]
  idx <- y == 0
  mean(y[idx] == yhat[idx])
}

BackErrorRate <- function(pred, yhat, y, at) {
  pred.order <- order(pred, decreasing = TRUE)
  y <- y[pred.order[1:at]]
  yhat <- yhat[pred.order[1:at]]
  mean(y != yhat)
}

BackPPV <- function(pred, yhat, y, at) {
  pred.order <- order(pred, decreasing = TRUE)
  y <- y[pred.order[1:at]]
  yhat <- yhat[pred.order[1:at]]
  idx <- yhat == 1
  mean(y[idx] == yhat[idx])
}

BackFMeasure <- function(pred, yhat, y, at, b) {
  prec <- BackPPV(pred, yhat, y, at)
  rec <- BackSensitivity(pred, yhat, y, at)
  2 * (prec * rec / (prec + rec))
}

BackAUC <- function(pred, yhat, y, at, b) {
  pred.order <- order(pred, decreasing = TRUE)
  y <- y[pred.order[1:at]]
  pred <- pred[pred.order[1:at]]
  as.numeric(pROC::auc(y, pred))
}