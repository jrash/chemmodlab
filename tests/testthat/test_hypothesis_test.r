context("Hypothesis Tests")

test_that("PerfCurveTest confidence intervals match those computed for Figure 5 in the paper", {
  
  skip_on_cran()
  
  load("test_hypothesis_test.rdata")
  
  # Simulate Case 1 scores from paper
  r <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
  r <- r[order(r)]/150000
  
  myseed <- 569
  set.seed(myseed)
  m <- 150000
  pi.0.true <- 1/501
  X <- rbinom(n=m, size=1, prob=pi.0.true)
  s.vec.ls <- list(c(1, 1), c(1, 1))
  params <- list(c(sqrt(2) * .8, 0), c(sqrt(2) * .6, 0))
  dist <- "binorm"
  rho <- .9
  # Parameters for bivariate normal distribution
  s1 <- s.vec.ls[[1]][1] ; s2 <- s.vec.ls[[2]][1];
  sigma1 <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2)
  s1 <- s.vec.ls[[1]][2]; s2 <- s.vec.ls[[2]][2];
  sigma2 <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2)# Covariance matrix
  bvn1 <- MASS::mvrnorm(m, mu = c(params[[1]][1], params[[2]][1]), Sigma = sigma1)*X + MASS::mvrnorm(m, mu = c(params[[1]][2], params[[2]][2]), Sigma = sigma2)*(1-X)
  S1 <- bvn1[, 1]
  S2 <- bvn1[, 2]
  
  method.v <- c("AH", "AH", "binomial", "binomial", "JZ Ind", "JZ Ind", "mcnemar", "mcnemar", "binomial Ind", "binomial Ind")
  correction.v <- c("none", "plus2", "none", "plus2", "none", "plus2", "none", "plus2", "none", "plus2")
  
  # JZ Ind is failing for some reason
  for(i in c(1:10)) {
    
    # Test recall
    ls <- PerfCurveTest(S1, S2, X, r, metric = "k", method = method.v[i], 
                        correction = correction.v[i], conf.level = .95)
    act_ci <- ls[[2]]
    
    exp_ci <- k.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
    # Test precision
    ls <- PerfCurveTest(S1, S2, X, r, metric = "pi", method = method.v[i], 
                        correction = correction.v[i],  conf.level = .95)
    act_ci <- ls[[2]]
    
    exp_ci <- pi.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
  }
  
})