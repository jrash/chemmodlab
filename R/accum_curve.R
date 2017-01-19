
HitCurve <- function(phat.list, y, max.select = NA, phat.labels = NA,
                     title = "") {
  # Order the cases by decreasing phat values, and plot the expected
  # number and actual number of hits as cases are selected.
  # Cases with tied phat values are grouped together.
  # phat       - list of estimated probabilities of a hit for test cases from
  #              some models
  # y          - numeric 0/1 indicator for the same cases giving the
  #              actual class (1 = hit)
  # max.select - maximum number of cases to be selected
  #
  #  modified from code originally written by
  #  William J. Welch 2001-2002{

  if (missing(phat.labels))
    phat.labels <- names(phat.list)

  list.index <- 1
  palette(c("black", "grey", rainbow(20)[c(1, 6, 11, 16, 2, 12, 17, 3, 13, 4, 19,
                                           10, 15, 1, 6, 11, 16, 2, 12, 17, 3,
                                           13, 4, 19, 10, 15)]))
  # palette(c('black','grey','blue','cyan','red','green','tomato','purple',
  # 'magenta','orange','seagreen','brown','yellow','steelblue',
  # 'blue','cyan','red','green','tomato','purple','magenta','orange','seagreen',
  # 'brown','yellow','steelblue'))
  ch.list <- c(45, 3, 0, 2, 4, 5, 6, 8, 11, 13, 12, 14)
  layout(matrix(c(1, 2), 1, 2), c(3.5, 1))

  # change any non-zero responses to 1
  y <- as.integer(y != 0)

  # if max.select not specified use min(300,#compounds/4)
  if (missing(max.select))
    max.select <- min(300, (length(y)/4))

  # This function will modify graphical parameters
  # Reset old parameters upon exiting
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  while (list.index <= length(phat.labels)) {
    if (length(phat.labels) == 1)
      phat <- unlist(phat.list)
    else phat <- unlist(phat.list[list.index])

      # Unique phat's, sorted with largest first.
      uniq.phat <- rev(sort(unique(phat)))

      select <- vector("numeric", length(uniq.phat))
      nhits <- vector("numeric", length(uniq.phat))

      for (i in 1:length(uniq.phat)) {
        cases.sel <- (phat == uniq.phat[i])
        select[i] <- sum(cases.sel)
        nhits[i] <- sum(y[cases.sel])
        if (sum(select[1:i]) >= max.select || i == length(uniq.phat)) {
          i.max <- i
          break
        }
      }

      uniq.phat <- uniq.phat[1:i.max]
      nhits <- nhits[1:i.max]
      select <- select[1:i.max]

      # Plot ideal HitCurve & random expectation the first time through
      if (list.index == 1) {
        x.max <- min(sum(select), max.select)
        y.max <- sum(y)
        par(mar = c(3, 3, 2, 0.5), mgp = c(3, 0.5, 0))
        plot(1:x.max, cumsum(rev(sort(y)))[1:x.max], type = "l",
             xlim = c(0, x.max), ylim = c(0, y.max), xlab = "", ylab = "",
             pch = "-")
        mtext(text = "Number of compounds selected", side = 1, line = 1.5)
        mtext(text = "Number of actual hits", side = 2, line = 1.5)
        title(title, cex = 0.8)
        lines(1:x.max, cumsum(rep(sum(y)/length(y), x.max)), lty = 1, col = 2)
      }

      # Plot the cumulative number of hits.
      cum.select.last <- select[1]
      cum.nhits.last <- nhits[1]
      for (i in 2:i.max) {
        cum.select <- cum.select.last + select[i]
        cum.nhits <- cum.nhits.last + nhits[i]
        # plot points and connect with lines
        # if(cum.select <= max.select){
          points(cum.select, cum.nhits, pch = ch.list[list.index + 1], cex = 0.3,
                 col = (list.index + 2))
          # }
        # if(cum.select <= max.select){
          lines(c(cum.select.last, cum.select), c(cum.nhits.last, cum.nhits),
                lty = 1, col = (list.index + 2))
          # }
        cum.select.last <- cum.select
        cum.nhits.last <- cum.nhits
      }

      list.index <- list.index + 1
  }
  par(mar = c(3, 0, 2, 0))

  plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
  legend("topright", c("Ideal", phat.labels, "Random"), cex = .75,
         bty = "n", title = "Classification",
         pch = c(ch.list[1:(list.index)], ch.list[1]),
         col = c(1, 3:(list.index + 1), 2),
         text.col = c(1, 3:(list.index + 1), 2))
}


#--------------------------------------------------------------------------------

ContCurve <- function(yhat.list, y, max.select = NA, yhat.labels = NA, title = "",
                      curves.only = FALSE, start.col = 0) {
  if (missing(yhat.labels))
    yhat.labels <- names(yhat.list)

  # if max.select not specified use min(300,#compounds/4)
  if (missing(max.select))
    max.select <- min(300, (length(y)/4))

  list.index <- 1
  palette(c("black", "grey", rainbow(20)[c(1, 6, 11, 16, 2, 12, 17, 3, 13, 4, 19,
                                           10, 15, 1, 6, 11, 16, 2, 12, 17, 3, 13,
                                           4, 19, 10, 15)]))
  # palette(c('black','grey','blue','cyan','red','green','tomato','purple','magenta',
  # 'orange','seagreen','brown','yellow','steelblue',
  # 'blue','cyan','red','green','tomato','purple','magenta','orange','seagreen',
  # 'brown','yellow','steelblue'))
  x.max <- min(length(y), max.select)
  y.max <- max(cumsum(rev(sort(y)))[1:x.max])
  y.min <- min(cumsum(rev(sort(y)))[1:x.max], 0)

  # This function will modify graphical parameters
  # Reset old parameters upon exiting
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  if (curves.only) {
    layout(matrix(c(1, 2), 1, 2), c(3.5, 1))
    par(mfg = c(1, 1))
    par(mar = c(3, 3, 2, 0.5), mgp = c(3, 0.5, 0))
    plot((1:x.max), cumsum(rev(sort(y)))[1:x.max], type = "n", ylab = "", xlab = "",
         xlim = c(0, x.max), ylim = c(y.min, y.max), axes = FALSE)
  } else {
    layout(matrix(c(1, 2), 1, 2), c(3.5, 1))
    par(mar = c(3, 3, 2, 0.5), mgp = c(3, 0.5, 0))
    plot((1:x.max), cumsum(rev(sort(y)))[1:x.max], type = "l", lty = 1,
         xlim = c(0, x.max), ylim = c(y.min, y.max), ylab = "", xlab = "")
    mtext(text = "Number of compounds selected", side = 1, line = 1.5)
    mtext(text = "Cumulative response", side = 2, line = 1.5)
    lines(1:x.max, cumsum(rep(sum(y)/length(y), x.max)), lty = 1, col = 2)
    title(title, cex = 0.8)
  }

  while (list.index <= length(yhat.labels)) {
    if (length(yhat.labels) == 1){
      yhat <- unlist(yhat.list)
    } else {
      yhat <- unlist(yhat.list[list.index])
    }
    order <- order(yhat, -y, decreasing = TRUE)
    ploty <- cumsum(y[order])[1:x.max]
    lines((1:x.max), ploty, lty = (((list.index - 1)%%5) + 2), col = list.index +
            2 + start.col)
    list.index <- list.index + 1
  }

  if (list.index > 1) {
    lty.list <- c(1)
    for (i in 1:(list.index - 1)) lty.list <- c(lty.list, (((i - 1)%%5) + 2))
    lty.list <- c(lty.list, 1)
    par(mar = c(3, 0, 2, 0))
    plot(1, 1, type = "n", xlab = "", ylab = "", axes = FALSE)
    legend("bottomright", c("Ideal", yhat.labels, "Random"), cex = .75, bty = "n",
           title = "Continuous", lty = lty.list,
           col = c(1, (3 + start.col):(list.index + 1 + start.col), 2),
           text.col = c(1, (3 + start.col):(list.index + 1 + start.col), 2))
  }
}
