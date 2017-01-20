
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
