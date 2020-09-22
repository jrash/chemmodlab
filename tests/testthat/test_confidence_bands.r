context("Confidence Bands")

test_that("PerfCurveBands pointwise confidence intervals match those computed for Case 1 in the paper", {
  
  skip_on_cran()
  # TODO make sure this works on your windows machine
  skip_on_os('mac')
  skip_on_os('linux')

  load("test_confidence_bands_pw_type.rdata")

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
  correction.v <- c(F, T, F, T, F, T)

  for(i in 1:6) {
    
    # Test recall
    act_ci <- PerfCurveBands(S, X, r, metric = "k", method = method.v[i], type = "pointwise",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100, myseed = 567)
    
    exp_ci <- k.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
    # Test precision
    act_ci <- PerfCurveBands(S, X, r, metric = "pi", method = method.v[i], type = "pointwise",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100, myseed = 567)
    
    
    exp_ci <- pi.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
  }
  
})


test_that("PerfCurveBands confidence bands match those computed for Case 1 in the paper", {
  
  skip_on_cran()
  
  load("test_confidence_bands_band_type.rdata")
  
  # Simulate Case 1 scores from paper
  r <- c(2^seq(1, 13), 3^seq(1, 8), 10, 105, 300, 1500, 15000)
  r <- r[order(r)]/150000
  
  m <- 150000
  pi.0.true <- 1/501
  set.seed(111)
  X <- rbinom(n=m, size=1, prob=pi.0.true)
  mu.mix <- c(1.4, 0)
  sd.mix <- 1
  S <- rnorm(n=m, mu.mix[1], sd.mix)*X + rnorm(n=m, mu.mix[2], sd.mix)*(1-X)
  
  method.v <- c("sup-t", "sup-t", "theta-proj", "theta-proj")
  correction.v <- c(F, T, F, T, F, T)
  
  for(i in 1:4) {
    
    # Test recall
    act_ci <- PerfCurveBands(S, X, r, metric = "k", method = method.v[i], type = "band",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100,
                             mc.rep = 100000, myseed = 111)
    
    exp_ci <- k.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
    # Test precision
    act_ci <- PerfCurveBands(S, X, r, metric = "pi", method = method.v[i], type = "band",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100,
                             mc.rep = 100000, myseed = 111)
    
    
    exp_ci <- pi.ls[[i]]$interval.ls[[1]]
    
    expect_equal(act_ci, exp_ci, check.attribute = F)
    
  }
  
})
