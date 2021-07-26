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
#' c("EmProc", "binomial", "JZ ind", "mcnemar", "binomial ind"). Precision options are
#' c("EmProc", "binomial", "JZ ind", "stouffer", "binomial ind").
#' @param plus should plus correction be used or not?
#' @param alpha the significance level.
#' 
#' @export
PerfCurveTest <- function(S1, S2, X, r, metric = "rec", method = "EmProc",
                          plus = T, alpha = .05, h = NULL, seed = 111){
  
  set.seed(seed)
  m <- length(S1) #total sample size
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
  hits1 <- hits2 <- hits12 <- ntest12 <- vector(length = length(r))
  for(i in 1:length(r)) {
    t1 <- F1inv(1-r[i])
    t2 <- F2inv(1-r[i])
    hits1[i] <- sum(X*(S1 > t1))
    hits2[i] <- sum(X*(S2 > t2))
    hits12[i] <- sum(X*(S1 > t1 & S2 > t2))
    ntest12[i] <- sum(S1 > t1 & S2 > t2)
  }
  nact <- sum(X) #total number of actives
  ntest <- (m*r) #number of compounds tested
  
  # Uncorrected parameters
  r <- ntest/m
  pi.0 <- nact/m
  k1 <- (hits1)/nact
  k2 <- (hits2)/nact
  k12 <- (hits12)/nact
  pi1 <- pi.0/r*k1
  pi2 <- pi.0/r*k2
  r12 <- ntest12/m
  pi12 <- hits12/ntest12
  pi12 <- ifelse(is.na(pi12), 0, pi12)
  # pooled estimates
  pi <- (pi1 + pi2)/2
  k <- (k1 + k2)/2
  
  if(plus) {
    hits1.c <- hits1 + 1
    hits2.c <- hits2 + 1
    nact.c <- nact + 2 
    ntest.c <- ntest + 1
    m.c <- m + 2
  } else {
    hits1.c <- hits1
    hits2.c <- hits2
    nact.c <- nact 
    ntest.c <- ntest
    m.c <- m
  }
  
  # Corrected parameters
  r.c <- ntest.c/m.c
  pi.0.c <- nact.c/m.c
  k1.c <- (hits1.c)/nact.c
  k2.c <- (hits2.c)/nact.c
  k12.c <- (hits12)/nact.c
  pi1.c <- (hits1.c)/(ntest.c)
  pi2.c <- (hits2.c)/(ntest.c)
  r12.c <- ntest12/m.c
  pi12.c <- hits12/ntest12
  pi12.c <- ifelse(is.na(pi12.c), 0, pi12.c)
  
  
  Lam1.vec <- vector(length = length(k1))
  Lam2.vec <- vector(length = length(k1))
  for(j in seq_along(k1)){
    Lam1.vec[j] <- EstLambda(S1, X, t = F1inv(1-r[j]), idx = idx[j], h)
    Lam2.vec[j] <- EstLambda(S2, X, t = F2inv(1-r[j]), idx = idx[j], h)
  }
  
  CI.int <- matrix(ncol = 2, nrow = length(k1))
  se.p <- se.np <- p.val <- vector(length = length(k1))
  if(metric == "rec") {
    if(method %in% c("JZ ind", "EmProc")){
      for(j in seq_along(k1)) {
        
        Lam1 <- Lam1.vec[j]
        Lam2 <- Lam2.vec[j]
        
        # unpooled variances
        var1.k.np <- ((k1.c[j]*(1-k1.c[j]))/(m.c*pi.0.c))*(1-2*Lam1) +
          (Lam1^2*(1-r.c[j])*r.c[j])/(m.c*pi.0.c^2)
        # Check to see if var.k.np is negative due to m.cachine precision problem.c
        var1.k.np <- ifelse(var1.k.np < 0, 0, var1.k.np)
        var2.k.np <- ((k2.c[j]*(1-k2.c[j]))/(m.c*pi.0.c))*(1-2*Lam2) +
          (Lam2^2*(1-r.c[j])*r.c[j])/(m.c*pi.0.c^2)
        var2.k.np <- ifelse(var2.k.np < 0, 0, var2.k.np)
        cov.k.np <- (m.c^-1*pi.0.c^-2)*(pi.0.c*(k12.c[j]-k1.c[j]*k2.c[j])*(1-Lam1-Lam2) +
                                          (r12.c[j]-r.c[j]^2)*Lam1*Lam2)
        
        # pooled variances          
        var1.k.p <- ((k[j]*(1-k[j]))/(m*pi.0))*(1-2*Lam1) +
          (Lam1^2*(1-r[j])*r[j])/(m*pi.0^2)
        var1.k.p <- ifelse(var1.k.p < 0, 0, var1.k.p)
        var2.k.p <- ((k[j]*(1-k[j]))/(m*pi.0))*(1-2*Lam2) +
          (Lam2^2*(1-r[j])*r[j])/(m*pi.0^2)
        var2.k.p <- ifelse(var2.k.p < 0, 0, var2.k.p)
        cov.k.p <- (m^-1*pi.0^-2)*(pi.0*(k12[j]-k[j]*k[j])*(1-Lam1-Lam2) +
                                     (r12[j]-r[j]^2)*Lam1*Lam2)
        
        if(method == "JZ ind"){
          # assuming independence of k
          var.k.p <- var1.k.p + var2.k.p
          var.k.np <- var1.k.np + var2.k.np
        } else if(method == "EmProc") {
          var.k.p <- var1.k.p + var2.k.p - 2*cov.k.p
          var.k.np <- var1.k.np + var2.k.np - 2*cov.k.np
        }
        var.k.p <- ifelse(var.k.p < 0, 0, var.k.p)
        var.k.np <- ifelse(var.k.np < 0, 0, var.k.np)
        se.p[j] <- sqrt(var.k.p)
        se.np[j] <- sqrt(var.k.np)
        
        # Use pooled variance for tests, and unpooled for CIs
        CI.int[j, ] <- c( (k1[j]-k2[j]) - quant*se.np[j],
                          (k1[j]-k2[j]) + quant*se.np[j])
        zscore <- ( (k1[j]-k2[j]) )/se.p[j]
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      }
    } else if(method %in% c("binomial ind", "binomial")) {
      for(j in seq_along(k1.c)) {
        
        # unpooled variances
        var1.k.np <- (m.c*pi.0.c)^-1*(k1.c[j])*(1-k1.c[j])
        var2.k.np <- (m.c*pi.0.c)^-1*(k2.c[j])*(1-k2.c[j])
        cov.k.np <- (m.c*pi.0.c)^-1*(k12.c[j] - k1.c[j]*k2.c[j])
        
        # pooled variances
        var1.k.p <- (m*pi.0)^-1*(k[j])*(1-k[j])
        var2.k.p <- (m*pi.0)^-1*(k[j])*(1-k[j])
        cov.k.p <- (m*pi.0)^-1*(k12[j] - k[j]*k[j])
        
        if(method == "binomial") {
          var.k.p <- var1.k.p + var2.k.p - 2*cov.k.p
          var.k.np <- var1.k.np + var2.k.np - 2*cov.k.np
        } else if(method == "binomial ind") {
          var.k.p <- var1.k.p + var2.k.p
          var.k.np <- var1.k.np + var2.k.np
        }
        var.k.p <- ifelse(var.k.p < 0, 0, var.k.p)
        var.k.np <- ifelse(var.k.np < 0, 0, var.k.np)
        se.p[j] <- sqrt(var.k.p)
        se.np[j] <- sqrt(var.k.np)
        CI.int[j, ] <- c( (k1[j]-k2[j]) - quant*se.np[j],
                          (k1[j]-k2[j]) + quant*se.np[j])
        zscore <- ( (k1[j]-k2[j]) )/se.p[j]
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      } 
    } else if(method == "mcnemar") {
      for(j in seq_along(k1.c)) {
        BminusC <- hits1[j] - hits2[j]
        BplusC <- hits1[j] + hits2[j] - 2*hits12[j]
        BminusC.c <- hits1.c[j] - hits2.c[j]
        BplusC.c <- hits1.c[j] + hits2.c[j] - 2*hits12[j]
        zscore <- ifelse(BplusC>0, BminusC/sqrt(BplusC), 0)
        p.val[j] <- 2*pnorm(-abs(zscore))
        se.np[j] <- sqrt(BplusC.c - BminusC.c^2/nact.c) / nact.c
        # TODO use the actual difference as the center point for the CI?
        # No, use actual BP corrected McNemar for a fair comparison
        CI.int[j,] <- c( BminusC.c/nact.c - quant*se.np[j], BminusC.c/nact.c + quant*se.np[j])
      }
    }
    # difference estimates
    est <- (k1 - k2)
  } else if(metric == "prec") {
    # TODO implement pooled and unpooled variances here
    if(method %in% c("JZ ind", "EmProc")) {
      for(j in seq_along(pi1.c)) {
        Lam1 <- Lam1.vec[j]
        var1.pi <- (pi1.c[j]*(1-pi1.c[j]))/(m*r[j]) + (1-r[j])*(pi1.c[j]-Lam1)^2/(m*r[j])
        var1.pi <- ifelse(var1.pi < 0, 0, var1.pi)
        Lam2 <- Lam2.vec[j]
        var2.pi <- (pi2.c[j]*(1-pi2.c[j]))/(m*r[j]) + (1-r[j])*(pi2.c[j]-Lam2)^2/(m*r[j])
        var2.pi <- ifelse(var2.pi < 0, 0, var2.pi)
        cov.pi <- ((m^-1*r[j]^-2))*(r12[j]*pi12.c[j]*(1 - pi12.c[j]) + 
                                      (pi12.c[j]-Lam1)*(pi12.c[j]-Lam2)*(r12[j] - r[j]^2))
        
        if(method == "JZ ind") {
          var.pi <- var1.pi + var2.pi
        } else if(method == "EmProc") {
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
      for(j in seq_along(k1.c)) {
        var1.pi <- (pi.c[j]*(1-pi.c[j]))/(m*r[j])
        var2.pi <- (pi.c[j]*(1-pi.c[j]))/(m*r[j])
        cov.pi <- (r12[j]/(m*r[j]^2))*(pi12.c[j]*(1 - pi12.c[j]))
        
        if(method == "binomial") {
          var.pi <- var1.pi + var1.pi - 2*cov.pi
        } else if(method == "binomial ind") {
          var.pi <- var1.pi + var1.pi
        }
        
        # Check to see if var.k is negative due to machine precision problem
        var.pi <- ifelse(var.pi < 0, 0, var.pi)
        se <- sqrt(var.pi)
        CI.int[j, ] <- cbind((pi1.c[j] - pi2.c[j]) - quant*se, (pi1.c[j] - pi2.c[j]) + quant*se)
        zscore <- (pi1.c[j] - pi2.c[j])/se
        p.val[j] <- 2*pnorm(-abs(zscore))
        p.val[j] <- ifelse(is.na(p.val[j]), 1, p.val[j])
      }
    } else if(method == "stouffer") {
      # TODO check that stouffer implemented correctly after changes
      for(j in seq_along(pi1.c)) {
        Sorder1 <-  order(S1,decreasing=TRUE)
        Sorder2 <-  order(S2,decreasing=TRUE)
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
        if(correction == "plus") {
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
  list(diff_estimate = est, std_err = se.np, ci_interval = CI.int, p_value = p.val)
}