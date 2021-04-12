#' Perform a hypothesis test for the difference between two performance curves.
#' 
#' \code{PerfCurveTest} takes score vectors for two scoring algorithms and an activity vector.
#' A performance curve is created for the two scoring algorithms and hypothesis tests are performed
#' at the selected testing fractions. 
#' 
#' @param S1 a vector of scores for scoring algorithm 1.
#' @param S2 a vector of scores for scoring algorithm 2.
#' @param X a vector of activities.
#' @param r a vector of testing fractions.
#' @param metric the performance curve to use. Options are recall ("rec") and precision ("prec").
#' @param method the method to use. Recall options are 
#' c("AH", "binomial", "JZ ind", "mcnemar", "binomial ind"). Precision options are
#' c("AH", "binomial", "JZ ind", "stouffer", "binomial ind").
#' @param alpha the significance level.
#' 
#' @export
PerfCurveTest <- function(S1, S2, X, r, metric = "rec", method = "AH",
                          correction = "plus2", alpha = .05, h = NULL){
  
  m <- length(S1) #total sample size
  nact <- sum(X)  #total number of actives
  r.all <- (1:m)/m
  idx <- which(r.all %in% r)
  
  # Z Quantile for CIs
  quant <- qnorm(1-alpha/2)
  
  # Convert fractions in r to thresholds, for both sets of scores
  # Also accumulate the number of actives, for each score and jointly
  F1cdf <- ecdf(S1); yyobs <- sort(unique(S1))
  F1inv <- stepfun(x = F1cdf(yyobs), y = c(yyobs, max(yyobs)), right=TRUE, f=1)
  F2cdf <- ecdf(S2); yyobs <- sort(unique(S2))
  F2inv <- stepfun(x = F2cdf(yyobs), y = c(yyobs, max(yyobs)), right=TRUE, f=1)
  hits1 <- hits2 <- hits12 <- r12 <- vector(length = length(r))
  for(i in 1:length(r)) {
    t1 <- F1inv(1-r[i])
    t2 <- F2inv(1-r[i])
    hits1[i] <- sum(X*(S1 > t1))
    hits2[i] <- sum(X*(S2 > t2))
    hits12[i] <- sum(X*(S1 > t1 & S2 > t2))
    r12[i] <- mean(S1 > t1 & S2 > t2)
  }
  
  pi.0 <- mean(X) #fraction of actives
  if(correction == "plus2") {
    # Treat the additional hits as if they were not same compounds
    # In either scoring algorithm.  The reason for this is that treating the
    # the compounds as shared can inflate the covariance when the number of 
    # compounds tested is small, leading to an increased type I error rate.
    # Want the plus 2 correction to have the maximum increase of the variance.
    pi1 <- (hits1+2)/(m*r+4)
    k1 <- r/pi.0*pi1
    pi2 <- (hits2+2)/(m*r+4)
    k2 <- r/pi.0*pi2
    pi <- (pi1 + pi2)/2
    k <- (k1 + k2)/2
    pi12 <- (hits12)/(m*r+4)
    k12 <- r/pi.0*pi12
  } else if(correction == "none") {
    pi1 <- (hits1)/(m*r)
    k1 <- (hits1)/nact
    pi2 <- (hits2)/(m*r)
    k2 <- (hits2)/nact
    pi12 <- (hits12)/(m*r12)
    k12 <- (hits12)/nact
    pi <- (pi1 + pi2)/2
    k <- (k1 + k2)/2
  }
  pi12 <- ifelse(is.na(pi12), 0, pi12)
  
  Lam1.vec <- vector(length = length(k))
  Lam2.vec <- vector(length = length(k))
  for(j in seq_along(k)){
    Lam1.vec[j] <- EstLambda(S1, X, t = F1inv(1-r[j]), idx = idx[j], h)
    Lam2.vec[j] <- EstLambda(S2, X, t = F2inv(1-r[j]), idx = idx[j], h)
  }
  
  CI.int <- matrix(ncol = 2, nrow = length(k1))
  se <- p.val <- vector(length = length(pi1))
  if(metric == "rec") {
    if(method %in% c("JZ ind", "AH")){
      for(j in seq_along(k1)) {
        Lam1 <- Lam1.vec[j]
        var1.k <- ((k1[j]*(1-k1[j]))/(m*pi.0))*(1-2*Lam1) +
          (Lam1^2*(1-r[j])*r[j])/(m*pi.0^2)
        # Check to see if var.k is negative due to machine precision problem
        var1.k <- ifelse(var1.k < 0, 0, var1.k)
        Lam2 <- Lam2.vec[j]
        var2.k <- ((k2[j]*(1-k2[j]))/(m*pi.0))*(1-2*Lam2) +
          (Lam2^2*(1-r[j])*r[j])/(m*pi.0^2)
        # Check to see if var.k is negative due to machine precision problem
        var2.k <- ifelse(var2.k < 0, 0, var2.k)
        cov.k <- (m^-1*pi.0^-2)*(pi.0*(k12[j]-k1[j]*k2[j])*(1-Lam1-Lam2) +
                                   (r12[j]-r[j]^2)*Lam1*Lam2)
        
        # assuming independence of k
        if(method == "JZ ind"){
          var.k <- var1.k + var2.k
        } else if(method == "AH") {
          var.k <- var1.k + var2.k - 2*cov.k
        }
        # Check to see if var.k is negative due to machine precision problem
        var.k <- ifelse(var.k < 0, 0, var.k)
        se[j] <- sqrt(var.k)
        CI.int[j, ] <- c( (hits1-hits2)[j]/nact - quant*se[j],    #jho change!
                          (hits1-hits2)[j]/nact + quant*se[j])
        zscore <- ( (hits1-hits2)[j]/nact )/se[j]    #jho change!
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      }
    } else if(method %in% c("binomial ind", "binomial")) {
      for(j in seq_along(k1)) {
        var1.k <- (m*pi.0)^-1*(k[j])*(1-k[j])
        var2.k <- (m*pi.0)^-1*(k[j])*(1-k[j])
        cov.k <- (m*pi.0)^-1*(k12[j] - k[j]*k[j])
        if(method == "binomial") {
          var.k <- var1.k + var1.k - 2*cov.k
        } else if(method == "binomial ind") {
          var.k <- var1.k + var1.k
        }
        var.k <- ifelse(var.k < 0, 0, var.k)
        se[j] <- sqrt(var.k)
        CI.int[j, ] <- c( (hits1-hits2)[j]/nact - quant*se[j],    #jho change!
                          (hits1-hits2)[j]/nact + quant*se[j])
        zscore <- ( (hits1-hits2)[j]/nact )/se[j]    #jho change!
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      } 
    } else if(method == "mcnemar") {
      for(j in seq_along(k1)) {
        mat <- matrix(c(hits12[j], hits1[j] - hits12[j],  hits2[j] - hits12[j],
                        nact - (hits1[j] + hits2[j] - hits12[j])), ncol = 2)
        if(!all(c(mat)>=0)) stop(paste(toString(mat),"Negative matrix?", 
                                       eff_seeds[z], toString(params), metric, 
                                       sep = "|"))
        BminusC <- hits1[j] - hits2[j]
        BplusC <- hits1[j] + hits2[j] - 2*hits12[j]
        zscore <- ifelse(BplusC>0, BminusC/sqrt(BplusC), 0)
        p.val[j] <- 2*pnorm(-abs(zscore))
        se[j] <- sqrt(BplusC - BminusC^2/nact) / nact
        CI.int[j,] <- c( BminusC/nact - quant*se[j], BminusC/nact + quant*se[j])
      }
    }
    # difference estimates
    est <- (hits1-hits2)/nact
  } else if(metric == "prec") {
    Sorder1 <-  order(S1,decreasing=TRUE)
    S1.o <- S1[Sorder1]
    c1 <- S1.o[idx]
    hits1 <- cumsum(X[Sorder1])[idx]
    Sorder2 <-  order(S2,decreasing=TRUE)
    S2.o <- S2[Sorder2]
    c2 <- S2.o[idx]
    hits2 <- cumsum(X[Sorder2])[idx]
    hits12 <- vector(length = length(idx))
    r12 <- vector(length = length(r))
    for(i in 1:length(idx)) {
      r12[i] <-  mean(S1 > c1[i] & S2 > c2[i])
      hits12[i] <- sum(X*(S1 > c1[i] & S2 > c2[i]))
      hits1[i] <- sum(X*(S1 > c1[i]))
      hits2[i] <- sum(X*(S2 > c2[i]))
    }
    
    if(method %in% c("JZ ind", "AH")) {
      Sorder1.idx <- Sorder1[idx]
      Sorder2.idx <- Sorder2[idx]
      for(j in seq_along(pi1)) {
        Lam1 <- Lam1.vec[j]
        var1.pi <- (pi1[j]*(1-pi1[j]))/(m*r[j]) + (1-r[j])*(pi1[j]-Lam1)^2/(m*r[j])
        var1.pi <- ifelse(var1.pi < 0, 0, var1.pi)
        Lam2 <- Lam2.vec[j]
        var2.pi <- (pi2[j]*(1-pi2[j]))/(m*r[j]) + (1-r[j])*(pi2[j]-Lam2)^2/(m*r[j])
        var2.pi <- ifelse(var2.pi < 0, 0, var2.pi)
        cov.pi <- ((m^-1*r[j]^-2))*(r12[j]*pi12[j]*(1 - pi12[j]) + 
                                      (pi12[j]-Lam1)*(pi12[j]-Lam2)*(r12[j] - r[j]^2))
        
        if(method == "JZ ind") {
          var.pi <- var1.pi + var2.pi
        } else if(method == "AH") {
          var.pi <- var1.pi + var2.pi - 2*cov.pi
        }
        # Check to see if var.k is negative due to machine precision problem
        var.pi <- ifelse(var.pi < 0, 0, var.pi)
        se <- sqrt(var.pi)
        CI.int[j, ] <- cbind((pi1[j] - pi2[j]) - quant*se, (pi1[j] - pi2[j]) + quant*se)
        zscore <- (pi1[j] - pi2[j])/se
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      }
    } else if(method %in% c("binomial ind", "binomial")) {
      for(j in seq_along(k1)) {
        var1.pi <- (pi[j]*(1-pi[j]))/(m*r[j])
        var2.pi <- (pi[j]*(1-pi[j]))/(m*r[j])
        cov.pi <- (r12[j]/(m*r[j]^2))*(pi12[j]*(1 - pi12[j]))
        
        if(method == "binomial") {
          var.pi <- var1.pi + var1.pi - 2*cov.pi
        } else if(method == "binomial ind") {
          var.pi <- var1.pi + var1.pi
        }
        # Check to see if var.k is negative due to machine precision problem
        var.pi <- ifelse(var.pi < 0, 0, var.pi)
        se <- sqrt(var.pi)
        CI.int[j, ] <- cbind((pi1[j] - pi2[j]) - quant*se, (pi1[j] - pi2[j]) + quant*se)
        zscore <- (pi1[j] - pi2[j])/se
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      }
    } else if(method == "stouffer") {
      for(j in seq_along(pi1)) {
        overlap.idx <- Sorder1[which(Sorder1[1:idx[j]] %in% Sorder2[1:idx[j]])]
        S1.unique.idx <- Sorder1[which(!Sorder1[1:idx[j]] %in% Sorder2[1:idx[j]])]
        S2.unique.idx <- Sorder2[which(!Sorder2[1:idx[j]] %in% Sorder1[1:idx[j]])]
        n12 <- length(overlap.idx)
        a <- sum(X[overlap.idx])
        d <- n12 - a
        n1 <- m*r[j] - n12 
        e <- sum(X[S1.unique.idx])
        g <- sum(X[S2.unique.idx])
        p.1 <- e/n1
        p.2 <- g/n1
        p.bar <- (e + g)/(2*n1)
        if(correction == "plus2") {
          n1 <- n1 + 4
          e <- e + 2
          g <- g + 2
        }
        z1 <- (p.1 - p.2)/sqrt((2*p.bar*(1-p.bar))/n1)
        z1 <- ifelse(is.na(z1), 0, z1)
        z3 <- 0
        w <- (n1 + n1)/(2*n12 + n1 + n1)
        z5 <- (w*z1 + (1-w)*z3)/(sqrt(w^2 + (1-w)^2))
        p.val[j] <- 2*(1 - pnorm(abs(z5)))
      }
    }
    # difference estimates
    est <- pi1 - pi2
  }
  list(diff_estimate = est, std_err = se, ci_interval = CI.int, p_value = p.val)
}
