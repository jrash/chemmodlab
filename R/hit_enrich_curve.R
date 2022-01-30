# Function to add shading ------
add.alpha <- function(col, alpha=1){
  if (missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha = alpha))  
}

#' Plot hit enrichment curves from observed scores and activities
#' 
#' Plot hit enrichment curves based on scores from multiple algorithms. Actual
#' activities are required. Additionally plot the ideal hit enrichment curve
#' that would result under perfect scoring, and the hit enrichment curve that
#' would result under random scoring. Optionally, simultaneous confidence bands
#' may also be requested.
#' 
#' @param S.df Data frame where variables are numeric scores from different
#'   algorithms. Rows represent unique compounds.
#' @param labels Character vector of labels for the different algorithms in
#'   \code{S.df}. If missing, variable names in \code{S.df} will be used.
#' @param y Numeric vector of activity values. Activity values must be either 0
#'   (inactive/undesirable) or 1 (active/desirable); no other values are
#'   accepted. Compounds are assumed to in the same order as in \code{S.df}.
#' @param x.max Integer, the maximum number of tests allowed on the x axis.
#' @param log Logical. \code{TRUE} plots the x axis on the log scale.
#' @param title Character string
#' @param conf Logical. \code{TRUE} plots (simultaneous) confidence bands for
#'   all hit enrichment curves.
#' @param conf.level Numeric, confidence coefficient
#' @param method Character indicates the method used to obtain confidence bands.
#'   The default is \code{sup-t} but other options (not recommended) are
#'   "theta-proj" and "bonf".
#' @param plus Logical. \code{TRUE} uses plus-adjusted version of \code{method}.
#' @param band.frac Numeric vector of fractions tested to be used in obtaining
#'   confidence bands. Vector should be no longer than \code{y}, and should have
#'   at least 20 entries. Entries should be in (0,1]. It is recommended that
#'   entries be consistent with between 1 and \code{x.max} tests.
#' 
#' @details By default, \code{x.max} is \code{length(y)}, so that hit enrichment
#' curves are obtained for all observable fractions, i.e., fractions of
#' \code{(1:length(y))/length(y)}. By default, confidence bands are evaluated
#' based on a smaller grid of 40 fractions. This smaller grid is evenly spaced
#' on either the original grid of \code{(1:length(y))/length(y)}, or the log
#' scale of the original grid.
#' 
#' @export
HitEnrich <-
  function(S.df,
           labels = NULL,
           y,
           x.max = NULL,
           log = TRUE, 
           title = "",
           conf = FALSE,
           conf.level = 0.95,
           method = "sup-t",
           plus = TRUE,
           band.frac = NULL
  ) {
    if (is.data.frame(S.df)) {
      if (is.null(labels)) labels <- colnames(S.df)
      m <- ncol(S.df)
    } else if (is.vector(S.df)) {
      if (is.null(labels)) labels <- "Score"
      m <- 1
      S.df <- cbind(S.df,NA)
    } else stop("S.df must be a data frame or a vector.")
    
    if (any((y != 0 & (y != 1)))) stop("y must be a vector of 0 and 1 only.")
    n <- length(y)
    n.actives <- sum(y)
    
    x.max <- ifelse(is.null(x.max), n, min(x.max, n))
    frac <- (1:x.max) / n
    if (log) frac <- log(frac, base = 10)
    
    safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", 
                                 "#332288", "#AA4499", "#44AA99", "#999933", 
                                 "#882255", "#661100", "#6699CC", "#888888")
    palette(safe_colorblind_palette)
    ch.list <- c(45, 3, 0, 2, 4, 5, 6, 8, 11, 13, 12, 14)
    
    if (conf) {
      if (missing(band.frac)) {
        band.frac <- frac[seq(1, x.max, len = 40)]
        if (log) { # spread fractions equally on log scale
          xvals <- seq(frac[1], frac[x.max], len = 40)
          tested <- unique(ceiling((10 ^ xvals) * n))
          band.frac <- log(tested / n, base = 10)
        }
      } else {
        band.frac <- unique(ceiling(band.frac * n)) / n
        if (log)
          band.frac <- log(band.frac, base = 10)
      }
      band.frac.idx <-  which(frac %in% band.frac)
    }
    
    mstart <- 0
    mleft <- m
    while (mleft > 0) {
      if (mleft <= 8) {
        for (idx in (mstart + 1):(mstart + mleft)) {
          if (idx == (mstart + 1)) { # plot the ideal and random curves
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
          
          if (conf) { # plot simultaneous confidence bands ----
            cl.ls <-
              PerfCurveBands(S = S, X = y, r = band.frac.idx / n, metric = "rec",
                             type = "band", method = method, plus = plus,
                             conf.level = conf.level)
            ucl <- cl.ls$CI[, 2]
            # Truncate the CI int at the ideal curve, because
            # the parameter cannot go above
            #ucl <- ifelse(ucl > y.ideal[band.frac.idx], y.ideal[band.frac.idx], ucl)
            lcl <- cl.ls$CI[, 1]
            polygon(c(band.frac, rev(band.frac)),
                    c(ucl, rev(lcl)),
                    col = add.alpha(idx, .25),
                    border = NA)
          }
          
          legend("bottomright", c("Ideal", labels[(mstart + 1):(mstart + mleft)], "Random"),
                 cex = 0.6, bty = "n", title = "Scoring Method",
                 pch = ch.list[1], 
                 col = c("black",1:mleft, "grey"), 
                 text.col = c("black",1:mleft, "grey"))
        }
      } else {
        for (idx in (mstart + 1):(mstart + 8)) {
          if (idx == (mstart + 1)) { # plot the ideal and random curves
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
          
          if (conf) { # plot simultaneous confidence bands ----
            cl.ls <-
              PerfCurveBands(S = S, X = y, r = band.frac.idx / n, metric = "rec",
                             type = "band", method = method, plus = plus,
                             conf.level = conf.level)
            ucl <- cl.ls$CI[, 2]
            # Truncate the CI int at the ideal curve, because
            # the parameter cannot go above
            #ucl <- ifelse(ucl > y.ideal[band.frac.idx], y.ideal[band.frac.idx], ucl)
            lcl <- cl.ls$CI[, 1]
            polygon(c(band.frac, rev(band.frac)),
                    c(ucl, rev(lcl)),
                    col = add.alpha(idx, .25),
                    border = NA)
          }
          
          legend("bottomright", c("Ideal", labels[(mstart + 1):(mstart + 8)], "Random"),
                 cex = 0.6, bty = "n", title = "Scoring Method",
                 pch = ch.list[1], 
                 col = c("black",1:8, "grey"), 
                 text.col = c("black",1:8, "grey"))
        }
        mstart <- mstart + 8
      }
      mleft <- mleft - 8
    }
    
  }

#' Plot differences between hit enrichment curves 
#' 
#' Plot differences between hit enrichment curves based on scores from multiple
#' algorithms. Actual activities are required. Additionally plot simultaneous
#' confidence bands for these differences. Plots may be used to determine if one
#' algorithm is "better" than another algorithm.
#' 
#' @param S.df Data frame where variables are numeric scores from at least 2
#'   different algorithms. Rows represent unique compounds.
#' @param labels Character vector of labels for the different algorithms in
#'   \code{S.df}. If missing, variable names in \code{S.df} will be used.
#' @param y Numeric vector of activity values. Activity values must be either 0
#'   (inactive/undesirable) or 1 (active/desirable); no other values are
#'   accepted. Compounds are assumed to in the same order as in \code{S.df}.
#' @param x.max Integer, the maximum number of tests allowed on the x axis.
#' @param log Logical. \code{TRUE} plots the x axis on the log scale.
#' @param title Character string
#' @param conf.level Numeric, confidence coefficient
#' @param method Character indicates the method used to obtain confidence bands.
#'   The default is \code{sup-t} but other options (not recommended) are
#'   "theta-proj" and "bonf".
#' @param plus Logical. \code{TRUE} uses plus-adjusted version of \code{method}.
#' @param band.frac Numeric vector of fractions tested to be used in obtaining
#'   confidence bands. Vector should be no longer than \code{y}, and should have
#'   at least 20 entries. Entries should be in (0,1]. It is recommended that
#'   entries be consistent with between 1 and \code{x.max} tests.
#' @param yrange Numeric vector of length 2. The desired range for the y axis.
#' 
#' @details By default, \code{x.max} is \code{length(y)}, so that hit enrichment
#' curves are obtained for all observable fractions, i.e., fractions of
#' \code{(1:length(y))/length(y)}. By default, confidence bands are evaluated
#' based on a smaller grid of 40 fractions. This smaller grid is evenly spaced
#' on either the original grid of \code{(1:length(y))/length(y)}, or the log
#' scale of the original grid.
#' 
#' @export
HitEnrichDiff <-
  function(S.df,
           labels = NULL,
           y,
           x.max = NULL,
           log = TRUE, 
           title = "",
           conf.level = 0.95,
           method = "sup-t",
           plus = TRUE,
           band.frac = NULL,
           yrange = NULL
  ) {
    if (is.data.frame(S.df)) {
      if (is.null(labels)) labels <- colnames(S.df)
      m <- ncol(S.df)
    } else stop("S.df must be a data frame with at least 2 columns of scores.")
    n <- length(y)
    #    n.actives <- sum(y)
    x.max <- ifelse(is.null(x.max), n, min(x.max, n))
    frac <- (1:x.max) / n
    if (log) frac <- log(frac, base = 10)
    
    if (missing(band.frac)) {
      band.frac <- frac[seq(1, x.max, len = 40)]
      if (log) {# spread fractions equally on log scale
        xvals <- seq(frac[1], frac[x.max], len = 40)
        tested <- unique(ceiling((10 ^ xvals) * n))
        band.frac <- log(tested / n, base = 10)
      }
    } else {
      band.frac <- unique(ceiling(band.frac * n)) / n
      if (log) band.frac <- log(band.frac, base = 10)
    }
    band.frac.idx <-  which(frac %in% band.frac)
    nt <- length(band.frac)
    
    K <- choose(m,2)
    setOfResults <- NULL
    score1 <- score2 <- diff_estim <- clb <- cub <- NULL
    k <- 0
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        k <- k + 1
        temp <- PerfCurveTest( S1 = S.df[, i], S2 = S.df[, j], 
                               X = y, r = band.frac.idx / n, metric = "rec", 
                               type = "band", method = method , plus = plus, 
                               pool = F, alpha = 1 - conf.level)
        score1 <- rbind(score1, labels[i])
        score2 <- rbind(score2, labels[j])
        diff_estim <- rbind(diff_estim, temp$diff_estimate)
        clb <- rbind(clb, temp$ci_interval[,1])
        cub <- rbind(cub, temp$ci_interval[,2])
      }
    }
    
    if (missing(yrange)) yrange <- c( min(clb), max(cub))
    xrange <- c( band.frac[1], band.frac[nt])
    
    par(mar = c(3, 3, 2, 0.5), mgp = c(3, 0.5, 0))
    for (k in 1:K) {
      plot(band.frac, NULL, xlim = xrange, ylim = yrange, type = "n", xaxt = "n")
      if (log) {
        axis(1, at = (-12:0), labels = 10 ^ (-12:0))
        abline(v = (-12:0), col = "lightgray", lty = "dotted")
        grid(nx = NA, ny = NULL)
        axis(3, at = seq(-12,0,.25), labels = round(10 ^ seq(-12,0,.25) * n),
             cex.axis = 0.5, tcl = NA)
      } else {
        axis(1)
        grid(nx = NULL, ny = NULL)
        axis(3, at = band.frac[seq(1, nt, len = 15)], 
             labels = round(band.frac * n)[seq(1, nt, len = 15)], 
             cex.axis = 0.5, tcl = NA)
      }
      abline(h = 0, col = "grey", lwd = 2)
      lines(band.frac, diff_estim[k,], col = "black", lty = 1, lwd = 1)
      polygon( c(band.frac, rev(band.frac)), c(cub[k,], rev(clb[k,])),
               col = add.alpha("black", .25), border = NA)
      mtext(text = "Fraction Tested", side = 1, line = 1.5)
      mtext(text = "Hit Enrichment Difference", side = 2, line = 1.5)
      mtext(text = "Number Tested", side = 3, line = 1.05, cex = 0.5)
      title(title, cex = 0.8)
      legend("topright", paste(score1[k], "-", score2[k]), cex = 0.6)
    }
    
  }