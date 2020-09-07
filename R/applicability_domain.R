#' Compute applicability domain for a chemmodlab model
#' 
#' \code{ApplicabilityDomain} evaluates the applicability domain for a chemmodlab model
#' using a Hotteling T2 control chart.
#' 
#' @param traindata training data
#' @param testdata test data
#' @param pvalue significance level for control limit threshold
#' @param desname descriptor set name
#' 
#' @import MSQC
#' @import MASS
#' 
#' @export
ApplicabilityDomain <- function(traindata, testdata, pvalue = .01,
                                desname = NULL) {
  
  # TODO Need a way to filter the same columns in the testdata that were
  # filtered by ModelTrain
    
  p1 <- HotellingControlChart(traindata, type = "t2", phase = 1,
                              alpha = pvalue)
  if (!is.null(desname))
    main.title <- paste0(desname, ": Hotelling Control Chart")
  else
    main.title <- "Hotelling Control Chart"
  #Phase II
  Xmv <- p1$Xmv
  S <- p1$covariance
  colm <- nrow(traindata)
  p2 <- HotellingControlChart(testdata, type = "t2", Xmv = Xmv, S = S,
                            colm = colm, main = main.title, phase = 2,
                            alpha = pvalue)
  test.outliers <- rownames(testdata)[which(p2$t2 > p2$ucl)]
  
  return(list(test.outliers = test.outliers))
}

HotellingControlChart <- function(type = c("chi", "t2", "mewma", "mcusum", "mcusum2"), 
                                   x, Xmv, S, colm, alpha = 0.01,
                                   k = 0.5, h = 5.5, phase = 1, method = "sw", main = NA, ...) {
  type <- match.arg(type)
  p <- ncol(x)
  m <- nrow(x)
  if (class(x) == "matrix" || class(x) == "data.frame") 
    (x <- array(data.matrix(x), c(m, p, 1)))
  n <- dim(x)[3]
  if (!missing(Xmv)) 
    (phase <- 2)
  x.jk <- matrix(0, m, p)
  t2 <- matrix(0, m, 1)
  x.jk <- apply(x, 1:2, mean)
  if (missing(Xmv)) 
    (Xmv <- colMeans(x.jk))
  if (missing(S)) 
    (S <- MSQC::covariance(x, method = method))
  if (missing(colm)) 
    (colm <- nrow(x))
  if (type == "t2") {
    for (ii in 1:m) {
      t2[ii, 1] <- n * t(x.jk[ii, ] - Xmv) %*% MASS::ginv(S) %*% 
        (x.jk[ii, ] - Xmv)
    }
    ifelse(n == 1, ifelse(phase == 1, ucl <- ((colm - 1)^2)/colm * 
      qbeta(1 - alpha, p/2, ((colm - p - 1)/2)), ucl <- ((p * 
      (colm + 1) * (colm - 1))/((colm^2) - colm * p)) * 
      qf(1 - alpha, p, colm - p)), ifelse(phase == 1, 
      ucl <- (p * (colm - 1) * (n - 1))/(colm * n - colm - 
        p + 1) * qf(1 - alpha, p, colm * n - colm - 
        p + 1), ucl <- (p * (colm + 1) * (n - 1))/(colm * 
        n - colm - p + 1) * qf(1 - alpha, p, colm * 
        n - colm - p + 1)))
  }
  t3 <- which(t2 > ucl)
  if (phase == 2) {
    par(mar = c(4, 5, 3, 5))
    plot(t2, ylim = c(0, 1.1 * max(max(t2), ucl)), main = main, 
      xlab = "Sample", ylab = expression(T^2), type = "o", 
      las = 1)
    points(t3, t2[t3], col = "red")
    segments(0, ucl, m, ucl, col = "red")
    mtext(paste(" UCL=", round(ucl, 2)), side = 4, at = ucl, 
      las = 2)
  }
  outList = list(main, ucl = round(ucl, 2), t2 = round(t2, 
    2), Xmv = round(Xmv, 2), covariance = signif(S, 2))
  return(outList)
}
