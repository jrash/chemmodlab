
Accumulation <- function(prob, y, at) {

  # Unique prob's, sorted with largest first.  if there are ties,
  # length(uniq.prob)<length(prob)
  uniq.prob <- sort(unique(prob), decreasing = TRUE)

  # select contains the number of cases at each unique prob nhits contains the
  # number of hits at each unique prob accpts contains the partial accumulation for
  # each prob (not unique)
  select <- vector("numeric", length(uniq.prob))
  nhits <- vector("numeric", length(uniq.prob))
  accpts <- vector("numeric", length(prob))

  accpts.index <- 1
  for (i in 1:length(uniq.prob)) {
    cases.sel <- (prob == uniq.prob[i])
    select[i] <- sum(cases.sel)
    nhits[i] <- sum(y[cases.sel])
    accpts[accpts.index:(accpts.index + select[i] - 1)] <- nhits[i]/select[i]
    accpts.index <- accpts.index + select[i]
  }

  # return the accumulation for every point, or those specified by 'at'
  if (missing(at))
    return(cumsum(accpts)) else return(cumsum(accpts)[at])
}


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
