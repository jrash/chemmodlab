# Returns Lambda, or P(+|S=t)
EstLambda = function(S, X, t, m){ 
  #h is the bandwidth, t is called c in the JZ paper
  #m is the sample size, S and X are the vectors of scores and labels
  h <- m^{-1/3}
  ind.win <- (S < t + h) & (S > t - h)
  exp.X.and.I <- sum(X*ind.win)/m
  exp.I <- sum(ind.win)/m
  exp.X.and.I/exp.I
}


# Bootstrap percentile CI
BootCI = function(X, S, m, pi.0, boot.rep, metric, plus2, r, myseed=111){
  storeout=matrix(NA, nrow=m, ncol=1+boot.rep)
  r.all <- (1:m)/m
  idx <- which(r.all %in% r)
  storeout <- storeout[idx, ]
  storeout[,1] <- r
  set.seed(myseed)
  for(i in 1:boot.rep){
    boot.samp <- sample(1:m, m, replace = T)
    X.star = X[boot.samp]
    S.star = S[boot.samp]
    Sorder <-  order(S.star,decreasing=TRUE)
    hits <- cumsum(X.star[Sorder])[idx] 
    pi.0 <- mean(X.star)
    pi <- (hits)/(m*r)
    k <- (hits)/(sum(X.star))
    if (plus2) {
      # using random plus2 correction
      plus2.yes <- rbinom(1, 4, .5)
      pi <- (hits+plus2.yes)/(m*r+4)
      k <- r/pi.0*pi
    }
    if(metric == "rec") {
      storeout[,1+i] = k
    } else if(metric == "prec") {
      storeout[,1+i] = pi
    } else if(metric == "lift") {
      storeout[,1+i] = k/r
    }
    # print(i/boot.rep)
  }
  A = apply(storeout[,-1],MARGIN=1,FUN=function(x){quantile(x, probs = c(.025, .975))})
  # Returns 95% percentile bootstrap intervals at quantiles
  A=t(A)
}


#' Construct a confidence band for a recall or precision curve
#' 
#' \code{PerfCurveBands} takes a pair of score and activity vectors as input.
#' A performance curve and confidence band is created for the selected testing fractions.
#' 
#' @param S a vector of scores.
#' @param X a vector of activities.
#' @param r a vector of testing fractions.
#' @param metric the performance curve to use. Options are recall ("rec") and precision ("prec").
#' @param type specifies whether a point-wise confidence interval 
#' ("pointwise") or a confidence band ("band") should be constructed.
#' @param method the method to use. Point-wise confidence interval options
#' are "binomial", "JZ", "bootstrap". Confidence band options are "sup-t", "theta-proj".
#' @param plus2 should plus2 correction be used or not?
#' @param conf.level the confidence level for the bands.
#' @param boot.rep the number of replicates to use for the bootstrap method.
#' @param mc.rep the number of Monte Carlo replicates to use for the sup-t method.
#' @param myseed the random seed.
#' 
#' 
#' @import MASS
#' @import signal
#' 
#' @export
#' 
PerfCurveBands <- function(S, X, r, metric = "rec", type = "band", method = "sup-t",
                           plus2 = T, conf.level = .95, boot.rep = 100,
                           mc.rep = 100000, myseed = 111){
  
  # Some error handeling
  
  if (!(method %in% c("binomial", "JZ", "bootstrap", "sup-t", "theta-proj"))) {
    stop("'method' should be a string specifiying a performance curve method in chemmodlab.
         Point-wise confidence interval options are 'binomial', 'JZ', 'bootstrap'. Confidence band options
         are 'sup-t', 'theta-proj'.")
  }
  
  if (!(metric %in% c("rec", "prec"))) {
    stop("'metric' should be a string specifiying a performance curve metric in chemmodlab.
         Point-wise confidence interval options are 'binomial', 'JZ', 'bootstrap'. Options
         are recall ('rec') and precision ('prec').")
  }
  
  if (!(type %in% c("band", "pointwise"))) {
    stop("'type' should be a string specifiying a performance curve type in chemmodlab.
         Point-wise confidence interval options are 'binomial', 'JZ', 'bootstrap'. Options
         are point-wise confidence interval ('pointwise') or a confidence band ('band').")
  }
  
  set.seed(myseed)
  alpha <- 1-conf.level
  
  # Compute indices of the testing fractions
  m <- length(S)
  r.all <- (1:m)/m
  # TODO This breaks when r contains zeros.  Need some error handling
  idx <- which(r.all %in% r)
  
  # Computing the performance measures
  Sorder <- order(S,decreasing=TRUE)
  Sorder.idx <- Sorder[idx]
  hits <- cumsum(X[Sorder])[idx]
  pi.0 <- mean(X)
  pi <- (hits)/(m*r)
  k <- (hits)/(sum(X))
  
  # Plus 2 correction
  if (plus2) {
    pi.c <- (hits+2)/(m*r+4)
    k.c <- r/pi.0*pi.c
  } else {
    pi.c <- pi
    k.c <- k
  }
  
  CI.int <- matrix(ncol = 2, nrow = length(k))
  if(type == "pointwise") {
    
    # Z Quantile for CIs
    quant <- qnorm(1-alpha/2)
    
    if(metric == "rec") {
      if(method == "JZ"){
        for(j in seq_along(k)) {
          Lam <- EstLambda(S[Sorder.idx], X[Sorder.idx], m, t = S[Sorder.idx][j])
          Lam <- signal::sgolayfilt(X[Sorder], p = 0, n = 11)[idx[j]]
          var.k <- ((k.c[j]*(1-k.c[j]))/(m*pi.0))*(1-2*Lam) + (Lam^2*(1-r[j])*r[j])/(m*pi.0^2)
          # Check to see if var.k is negative due to machine precision problem
          var.k <- ifelse(var.k < 0, 0, var.k)
          sd.k <- sqrt(var.k)
          CI.int[j, ] <- c(k[j] - quant*sd.k, k[j] + quant*sd.k)
        }
      } else if(method == "bootstrap") {
        # bootstrap quantiles
        CI.int <- BootCI(X, S, m, pi.0, boot.rep, metric = "rec", plus2, r, myseed=myseed)
      } else if(method == "binomial") {
        for(j in seq_along(k)) {
          var.k <- ((m*pi.0)^-1)*k.c[j]*(1-k.c[j])
          var.k <- ifelse(var.k < 0, 0, var.k)
          sd.k <- sqrt(var.k)
          CI.int[j, ]  <- cbind(k[j] - quant*sd.k, k[j] + quant*sd.k)
        }
      }
    } else if(metric == "prec") {
      if(method == "JZ") {
        for(j in seq_along(k)) {
          Lam <- EstLambda(S, X, m, t = S[Sorder.idx][j])
          Lam <- signal::sgolayfilt(X[Sorder], p = 0, n = 11)[idx[j]]
          var.pi <- (pi.c[j]*(1-pi.c[j]))/(m*r[j]) + (1-r[j])*(pi.c[j]-Lam)^2/(m*r[j])
          # Check to see if var.pi is negative due to machine precision issues
          var.pi <- ifelse(var.pi < 0, 0, var.pi)
          sd.pi <- sqrt(var.pi)
          CI.int[j, ] <- cbind(pi[j] - quant*sd.pi, pi[j] + quant*sd.pi)
        }
      } else if(method == "bootstrap") {
        # bootstrap quantiles
        CI.int <- BootCI(X, S, m, pi.0, boot.rep, metric = "prec", plus2, r, myseed=myseed)
      } else if(method == "binomial") {
        for(j in seq_along(pi)) {
          var.pi <- ((m*r[j])^-1)*pi.c[j]*(1-pi.c[j])
          # Check to see if var.pi is negative due to machine precision issues
          var.pi <- ifelse(var.pi < 0, 0, var.pi)
          sd.pi <- sqrt(var.pi)
          CI.int[j, ] <- cbind(pi[j] - quant*sd.pi, pi[j] + quant*sd.pi)
        }
      }
    }
  } else if(type == "band") {
      if(metric == "rec") {
        if(method == "sup-t") {
          cor.C <- matrix(NA, ncol = length(k), nrow = length(k))
          for(f in seq_along(k)) {
            for(e in 1:f) {
              Lam1 <- EstLambda(S, X, m, t = S[Sorder.idx][e])
              Lam1 <- signal::sgolayfilt(X[Sorder], p = 0, n = 11)[idx[e]]
              var.k1 <- ((k.c[e]*(1-k.c[e]))/(m*pi.0))*(1-2*Lam1) 
                        + (Lam1^2*(1-r[e])*r[e])/(m*pi.0^2)
              var.k1 <- ifelse(var.k1 < 0, 0, var.k1)
              Lam2 <- EstLambda(S, X, m, t = S[Sorder.idx][f])
              Lam2 <- signal::sgolayfilt(X[Sorder], p = 0, n = 11)[idx[f]]
              var.k2 <- ((k.c[f]*(1-k.c[f]))/(m*pi.0))*(1-2*Lam2) 
                        + (Lam2^2*(1-r[f])*r[f])/(m*pi.0^2)
              var.k2 <- ifelse(var.k2 < 0, 0, var.k2)
              cov.k <- (m^-1*pi.0^-2)*(pi.0*(k.c[e]-k.c[e]*k.c[f])*(1-Lam1-Lam2) 
                        + (r[e]-r[e]*r[f])*Lam1*Lam2)
              cov.k <- ifelse(cov.k < 0, 0, cov.k)
              cor.k <- cov.k/(sqrt(var.k1)*sqrt(var.k2))
              cor.k <- ifelse(var.k1 == 0 | var.k2 == 0, 0, cor.k)
              cor.C[e, f] <- cor.k
              cor.C[f, e] <- cor.k
            }
          }
          mc.samples <- MASS::mvrnorm(n = mc.rep, rep(0, length = length(k)), cor.C, tol = 1)
          max.q <- vector(length = mc.rep)
          for(j in 1:mc.rep) {
            max.q[j] <- max(abs(mc.samples[j, ]))
          }
          quant <- quantile(max.q, probs = 1-alpha/2)
        } else if(method == "theta-proj") {
          quant <- sqrt(qchisq(1-alpha, length(k)))
        }
        for(j in seq_along(k)) {
          Lam <- EstLambda(S, X, m, t = S[Sorder.idx][j])
          Lam <- signal::sgolayfilt(X[Sorder], p = 0, n = 11)[idx[j]]
          var.k <- ((k.c[j]*(1-k.c[j]))/(m*pi.0))*(1-2*Lam) + (Lam^2*(1-r[j])*r[j])/(m*pi.0^2)
          # Check to see if var.k is negative due to machine precision problem
          var.k <- ifelse(var.k < 0, 0, var.k)
          sd.k <- sqrt(var.k)
          CI.int[j, ] <- c(k[j] - quant*sd.k, k[j] + quant*sd.k)
        }
      } else if(metric == "prec") {
        if(method == "sup-t") {
          cor.C <- matrix(NA, ncol = length(k), nrow = length(k))
          for(f in seq_along(pi)) {
            for(e in 1:f) {
              Lam1 <- EstLambda(S, X, m, t = S[Sorder.idx][e])
              var.pi1 <- (pi.c[e]*(1-pi.c[e]))/(m*r[e]) + (1-r[e])*(pi.c[e]-Lam1)^2/(m*r[e])
              var.pi1 <- ifelse(var.pi1 < 0, 0, var.pi1)
              Lam2 <- EstLambda(S, X, m, t = S[Sorder.idx][f])
              var.pi2 <- (pi.c[f]*(1-pi.c[f]))/(m*r[f]) + (1-r[f])*(pi.c[f]-Lam2)^2/(m*r[f])
              var.pi2 <- ifelse(var.pi2 < 0, 0, var.pi2)
              cov.pi <- ((m*r[e]*r[f])^{-1})*(r[e]*pi.c[e]*(1 - pi.c[e])
                        + (pi.c[e]-Lam1)*(pi.c[e]-Lam2)*(r[e] - r[e]*r[f]))
              cov.pi <- ifelse(cov.pi < 0, 0, cov.pi)
              cor.pi <- cov.pi/(sqrt(var.pi1)*sqrt(var.pi2))
              cor.pi <- ifelse(var.pi1 == 0 | var.pi2 == 0, 0, cor.pi)
              cor.C[e, f] <- cor.pi
              cor.C[f, e] <- cor.pi
            }
          }
          mc.samples <- MASS::mvrnorm(n = mc.rep, rep(0, length = length(k)), cor.C, tol = 1)
          max.q <- vector(length = mc.rep)
          for(j in 1:mc.rep) {
            max.q[j] <- max(abs(mc.samples[j, ]))
          }
          quant <- quantile(max.q, probs = 1-alpha/2)
        } else if(method == "theta-proj") {
          quant <- sqrt(qchisq(1-alpha, length(pi)))
        }
        for(j in seq_along(pi)) {
          Lam <- EstLambda(S, X, m, t = S[Sorder.idx][j])
          var.pi <- (pi.c[j]*(1-pi.c[j]))/(m*r[j]) + (1-r[j])*(pi.c[j]-Lam)^2/(m*r[j])
          # Check to see if var.k is negative due to machine precision problem
          var.pi <- ifelse(var.pi < 0, 0, var.pi)
          sd.pi <- sqrt(var.pi)
          CI.int[j, ] <- cbind(pi[j] - quant*sd.pi, pi[j] + quant*sd.pi)
        }
      }
  }
  
  list(CI = CI.int, rec = k, prec = pi)
  
}
