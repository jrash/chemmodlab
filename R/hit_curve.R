#' Plot hit enrichment curves from observed scores and activities
#' 
#' Plot hit enrichment curves based on scores from multiple algorithms. Actual
#' activities are required. Additionally plot the ideal hit enrichment curve
#' that would result under perfect scoring, and the hit enrichment curve that
#' would result under random scoring. Optionally, simultaneous confidence bands
#' may also be requested.
#' 
#' @param S.df Data frame where variables are scores from different algorithms.
#'   Rows represent unique compounds.
#' @param labels Character vector of labels for the different algorithms in
#'   \code{S.df}. Variable names in \code{S.df} will be used by default.
#' @param y Numeric vector of non-negative activity values. Compounds are
#'   assumed to in the same order as in \code{S.df}. The value 0 is assumed to
#'   be undesirable. If \code{y.bin} is \code{TRUE}, then \code{y} will be
#'   converted using \code{y <- as.integer(y != 0)}.
#' @param y.bin Logical. \code{TRUE} converts \code{y} to a vector or 0 and 1.
#' @param x.max Integer, the maximum number of tests allowed on the x axis.
#' @param log Logical. \code{TRUE} plots the x axis on the log scale.
#' @param title Character string
HitCurve <-
  function(S.df,
           labels = NULL,
           y,
           y.bin = TRUE,
           x.max = NULL,
           log = TRUE, 
           title = ""
  ) {
    if (is.null(labels)) labels <- colnames(S.df)
    m <- ncol(S.df)
    if (y.bin) y <- as.integer(y != 0)
    n <- length(y)
    n.actives <- sum(y)
    # x.max <- ifelse(is.null(x.max), min(300, n / 4), min(x.max, n))
    x.max <- ifelse(is.null(x.max), n, min(x.max, n))
    frac <- (1:x.max) / n
    if (log) frac <- log(frac, base = 10)
    safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", 
                                 "#332288", "#AA4499", "#44AA99", "#999933", 
                                 "#882255", "#661100", "#6699CC", "#888888")
    palette(safe_colorblind_palette)
    ch.list <- c(45, 3, 0, 2, 4, 5, 6, 8, 11, 13, 12, 14)
    if (conf) {
      tested <- unique(ceiling(seq(frange[1], frange[2], len = fnum) * n))
      ntested <- length(tested)
      if (missing(band.frac)) {
        band.frac <- frac[seq(1, x.max, len = 40)]
      } else {
        band.frac <- unique(ceiling(band.frac * n))
        if (log) band.frac <- log(band.frac, base = 10)
      }
      band.frac.idx <-  which(frac %in% band.frac)
    }
    
    for (idx in 1:m) {
      if (idx == 1) { # plot the ideal and random curves
        y.ideal <- cumsum(rev(sort(y)))[1:x.max] / n.actives
        y.random <- cumsum(rep(n.actives / n, x.max)) / n.actives
        y.max <- max(y.ideal)
        par(mar = c(3, 3, 2, 0.5), mgp = c(3, 0.5, 0))
        
        plot(frac, y.ideal, col = "black", 	# ideal curve
             type = "l", xlab = "", ylab = "", pch = "-", xaxt = 'n')
        if (log) {
          axis(1, at = (-12:0), labels = 10 ^ (-12:0))
          abline(v = (-12:0), col = "lightgray", lty = "dotted")
          grid(nx = NA, ny = NULL)
          axis(3, at = seq(-12,0,.25), labels = round(10 ^ seq(-12,0,.25) * n),
               cex.axis = 0.5, tcl = NA)
        } else {
          axis(1)
          grid(nx = NULL, ny = NULL)
          axis(3, at = frac[seq(1, x.max, len = 15)], 
               labels = round(frac * n)[seq(1, x.max, len = 15)], 
               cex.axis = 0.5, tcl = NA)
        }
        mtext(text = "Fraction Tested", side = 1, line = 1.5)
        mtext(text = "Hit Enrichment", side = 2, line = 1.5)
        mtext(text = "Number Tested", side = 3, line = 1.05, cex = 0.5)
        title(title, cex = 0.8)
        lines(frac, y.random, lty = 1, col = "grey")	# random curve
      }
      
      # plot hit enrichment curves from algorithm scores ----
      S <- S.df[,idx]
      Fcdf <- ecdf(S); yyobs <- sort(unique(S))
      Finv <- stepfun(x = Fcdf(yyobs), y = c(yyobs, max(yyobs)), right = TRUE, 
                      f = 1)
      Fplus <- ecdf(S[y == 1])
      hits <- vector(length = x.max)
      for (i in 1:x.max) {
        t <- Finv(1 - i/n)
        hits[i] <- 1 - Fplus(t)
      }
      lines( frac, hits, lty = 1, col = idx, lwd = 2)
      
    }
    legend("bottomright", c("Ideal", labels, "Random"),
           cex = 0.6, bty = "n", title = "Scoring Method",
           pch = ch.list[1], 
           col = c("black",1:m, "grey"), 
           text.col = c("black",1:m, "grey"))
  }