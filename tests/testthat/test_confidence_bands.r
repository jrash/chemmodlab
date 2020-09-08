context("Confidence Bands - Pointwise")

test_that("PerfCurveBands confidence bands match those computed for Case 1 in the paper", {
  
  skip_on_cran()

  load("test_confidence_bands.rdata")

  # Simulate Case 1 scores from paper
  r <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
  r <- r[order(r)]/150000

  m <- 150000
  pi.0.true <- 1/501
  set.seed(567)
  X <- rbinom(n=m, size=1, prob=pi.0.true)
  mu.mix <- c(1.4, 0)
  sd.mix <- 1
  S <- rnorm(n=m, mu.mix[1], sd.mix)*X + rnorm(n=m, mu.mix[2], sd.mix)*(1-X)

  method.v <- c("binomial", "binomial", "JZ", "JZ", "bootstrap", "bootstrap")
  correction.v <- c("none", "JZ", "none", "JZ", "none", "JZ")

  for(i in 1:6) {
    
    # Test recall
    act_ci <- PerfCurveBands(S, X, r, metric = "k", method = method.v[i], correction = correction.v[i],
                             conf.level = .95, boot.rep = 100, myseed = 567)
    
    exp_ci <- k.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
    # Test precision
    act_ci <- PerfCurveBands(S, X, r, metric = "pi", method = method.v[i], correction = correction.v[i],
                             conf.level = .95, boot.rep = 100, myseed = 567)
    
    exp_ci <- pi.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
  }
  
})