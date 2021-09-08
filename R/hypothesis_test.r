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
#' @param plus should plus correction be applied to the confidence intervals?
#' @param pool use pooling for hypothesis tests? Only relevant to "EmProc".
#' @param alpha the significance level.
#' 
#' @export
PerfCurveTest <- function(S1, S2, X, r, metric = "rec", method = "EmProc", type = "pointwise", 
                          plus = T, pool = F, alpha = .05, h = NULL, seed = 111, mc.rep = 100000){
  
  # TODO need some more error checking here
  if(pool & method != "EmProc") warning("Pooling only relevant to EmProc and will not be applied")
  
  if (!(method %in% c("EmProc", "binomial", "JZ ind", "mcnemar", "binomial ind", "sup-t"))) {
    stop("'method' should be a string specifiying a performance curve method in chemmodlab.
         Point-wise confidence interval options are 'EmProc', 'binomial', 'JZ ind', 'mcnemar', 'binomial ind'. Confidence band options
         are 'sup-t'.")
  }
  
  set.seed(seed)
  m <- length(S1) #total sample size
  r.all <- (1:m)/m
  idx <- which(r.all %in% r)
  
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
    if(type == "pointwise") {
      quant <- qnorm(1-alpha/2)
      
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
          CI.int[j, ] <- c( (k1.c[j]-k2.c[j]) - quant*se.np[j],
                            (k1.c[j]-k2.c[j]) + quant*se.np[j])
          zscore <- ( (k1[j]-k2[j]) )/ifelse(pool, se.p[j], se.np[j])
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
          CI.int[j, ] <- c( (k1.c[j]-k2.c[j]) - quant*se.np[j],
                            (k1.c[j]-k2.c[j]) + quant*se.np[j])
          # For CorrBinom, unpooled is always used
          zscore <- ( (k1[j]-k2[j]) )/se.np[j]
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
          CI.int[j,] <- c( BminusC.c/nact.c - quant*se.np[j], BminusC.c/nact.c + quant*se.np[j])
        }
      }
    } else if (type == "band") {
      if(method == "sup-t") {
        cor.C <- matrix(NA, ncol = length(k1), nrow = length(k1))
        for(f in seq_along(k1)) {
          for(e in 1:f) {
            
            # Alg 1 unpooled variance (only using unpooled variances for bands)
            Lam11 <- Lam1.vec[e]; Lam21 <- Lam2.vec[e];
            var1.k.r1 <- ((k1.c[e]*(1-k1.c[e]))/(m.c*pi.0.c))*(1-2*Lam11) +
              (Lam11^2*(1-r.c[e])*r.c[e])/(m.c*pi.0.c^2)
            var1.k.r1 <- ifelse(var1.k.r1 < 0, 0, var1.k.r1)
            var2.k.r1 <- ((k2.c[e]*(1-k2.c[e]))/(m.c*pi.0.c))*(1-2*Lam21) +
              (Lam21^2*(1-r.c[e])*r.c[e])/(m.c*pi.0.c^2)
            var2.k.r1 <- ifelse(var2.k.r1 < 0, 0, var2.k.r1)
            cov.k.r1 <- (m.c^-1*pi.0.c^-2)*(pi.0.c*(k12.c[e]-k1.c[e]*k2.c[e])*(1-Lam11-Lam21) +
                                              (r12.c[e]-r.c[e]^2)*Lam11*Lam21)
            var.difk.r1 <- var1.k.r1 + var2.k.r1 - 2*cov.k.r1
            
            # Alg 2 unpooled variances
            Lam12 <- Lam1.vec[f]; Lam22 <- Lam2.vec[f];
            var1.k.r2 <- ((k1.c[f]*(1-k1.c[f]))/(m.c*pi.0.c))*(1-2*Lam12) +
              (Lam12^2*(1-r.c[f])*r.c[f])/(m.c*pi.0.c^2)
            var1.k.r2 <- ifelse(var1.k.r2 < 0, 0, var1.k.r2)
            var2.k.r2 <- ((k2.c[f]*(1-k2.c[f]))/(m.c*pi.0.c))*(1-2*Lam22) +
              (Lam22^2*(1-r.c[f])*r.c[f])/(m.c*pi.0.c^2)
            var2.k.r2 <- ifelse(var2.k.r2 < 0, 0, var2.k.r2)
            cov.k.r2 <- (m.c^-1*pi.0.c^-2)*(pi.0.c*(k12.c[f]-k1.c[f]*k2.c[f])*(1-Lam12-Lam22) +
                                              (r12.c[f]-r.c[f]^2)*Lam12*Lam22)
            var.difk.r2 <- var1.k.r2 + var2.k.r2 - 2*cov.k.r2
            
            # Covariance
            cov.k11.k12 <- ((m.c^-1*pi.0.c^-2)*(pi.0.c*(k1.c[e]-k1.c[e]*k1.c[f])*(1-Lam11-Lam12) +
                                                  (r.c[e]-r.c[e]*r.c[f])*Lam11*Lam12))
            cov.k21.k22 <- ((m.c^-1*pi.0.c^-2)*(pi.0.c*(k2.c[e]-k2.c[e]*k2.c[f])*(1-Lam21-Lam22) +
                                                  (r.c[e]-r.c[e]*r.c[f])*Lam21*Lam22))
            
            t1 <- F1inv(1-r[e])
            t2 <- F2inv(1-r[f])
            hits11.22 <- sum(X*(S1 > t1 & S2 > t2))
            ntest11.22 <- sum(S1 > t1 & S2 > t2)
            k11.22.c <- (hits11.22)/nact.c
            r11.22.c <- ntest11.22/m.c
            
            cov.k11.k22 <- (m.c^-1*pi.0.c^-2)*(pi.0.c*(k11.22.c-k1.c[e]*k2.c[f])*(1-Lam11-Lam22) +
                                                 (r11.22.c-r.c[e]*r.c[f])*Lam11*Lam22)
            
            t1 <- F1inv(1-r[f])
            t2 <- F2inv(1-r[e])
            hits12.21 <- sum(X*(S1 > t1 & S2 > t2))
            ntest12.21 <- sum(S1 > t1 & S2 > t2)
            k12.21.c <- (hits12.21)/nact.c
            r12.21.c <- ntest12.21/m.c
            
            cov.k12.k21 <- (m.c^-1*pi.0.c^-2)*(pi.0.c*(k12.21.c-k1.c[f]*k2.c[e])*(1-Lam21-Lam12) +
                                                 (r12.21.c-r.c[f]*r.c[e])*Lam21*Lam12)
            
            cov.difk <- cov.k11.k12 + cov.k21.k22 - cov.k11.k22 - cov.k12.k21
            
            if(e == f) {
              # If the covariance serves as a variance then it cant be non-negative
              # otherwise, it can be negative
              cov.difk <- ifelse(cov.difk < 0, 0, cov.difk)
            }
            cor.difk <- cov.difk/(sqrt(var.difk.r1)*sqrt(var.difk.r2))
            cor.difk <- ifelse(var.difk.r1 == 0 | var.difk.r2 == 0, ifelse(e == f, 1, 0), cor.difk)
            cor.C[e, f] <- cor.difk
            cor.C[f, e] <- cor.difk
          }
        }
        mc.samples <- MASS::mvrnorm(n = mc.rep, rep(0, length = length(k)), cor.C, tol = 1)
        max.q <- vector(length = mc.rep)
        for(j in 1:mc.rep) {
          max.q[j] <- max(abs(mc.samples[j, ]))
        }
        quant <- quantile(max.q, probs = 1-alpha)
      } else if(method == "theta-proj") {
        quant <- sqrt(qchisq(1-alpha, length(k1)))
      } else if(method == "bonf") {
        quant <- qnorm(1-alpha/(2*length(k1)))
      }
      for(j in seq_along(k)) {
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
        var.k.np <- var1.k.np + var2.k.np - 2*cov.k.np
        var.k.np <- ifelse(var.k.np < 0, 0, var.k.np)
        se.np[j] <- sqrt(var.k.np)
        CI.int[j, ] <- c( (k1.c[j]-k2.c[j]) - quant*se.np[j],
                          (k1.c[j]-k2.c[j]) + quant*se.np[j])
        p.val[j] <- NA
      }
    }
    # difference estimates
    est <- (k1 - k2)
  }
  list(diff_estimate = est, std_err = se.np, ci_interval = CI.int, p_value = p.val)
}