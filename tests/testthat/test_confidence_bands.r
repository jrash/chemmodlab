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

  r.all <- (1:m)/m
  idx <- which(r.all %in% r)
  k.ide <- (cumsum(rev(sort(X)))[idx]/sum(X))
  
  for(i in 1:6) {
    
    # Test recall
    act_ci <- PerfCurveBands(S, X, r, metric = "rec", method = method.v[i], type = "pointwise",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100, myseed = 567)
    
    exp_ci <- k.ls[[i]]$interval.ls[[1]]
    exp_ci[, 1] <- ifelse(exp_ci[, 1] < 0, 0, exp_ci[, 1])
    exp_ci[, 2] <- ifelse(exp_ci[, 2] > k.ide, k.ide, exp_ci[, 2])
    
    expect_equal(act_ci$CI, exp_ci, check.attribute = F)
    
    # Test precision
    act_ci <- PerfCurveBands(S, X, r, metric = "prec", method = method.v[i], type = "pointwise",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100, myseed = 567)
    
    
    exp_ci <- pi.ls[[i]]$interval.ls[[1]]
    exp_ci[, 1] <- ifelse(exp_ci[, 1] < 0, 0, exp_ci[, 1])
    exp_ci[, 2] <- ifelse(exp_ci[, 2] > 1, 1, exp_ci[, 2])
    
    expect_equal(act_ci$CI, exp_ci, check.attribute = F)
    
  }
  
})


test_that("PerfCurveBands confidence bands match those computed for Case 1 in the paper", {
  
  skip_on_cran()
  # TODO make sure this works on your windows machine
  skip_on_os('mac')
  skip_on_os('linux')
  
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
    act_ci <- PerfCurveBands(S, X, r, metric = "rec", method = method.v[i], type = "band",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100,
                             mc.rep = 100000, myseed = 111)
    
    exp_ci <- k.ls[[i]]$interval.ls[[1]]
    exp_ci[, 1] <- ifelse(exp_ci[, 1] < 0, 0, exp_ci[, 1])
    exp_ci[, 2] <- ifelse(exp_ci[, 2] > k.ide, k.ide, exp_ci[, 2])
    
    expect_equal(act_ci$CI, exp_ci, check.attribute = F)
    
    # Test precision
    act_ci <- PerfCurveBands(S, X, r, metric = "prec", method = method.v[i], type = "band",
                             plus2 = correction.v[i], conf.level = .95, boot.rep = 100,
                             mc.rep = 100000, myseed = 111)
    
    
    exp_ci <- pi.ls[[i]]$interval.ls[[1]]
    exp_ci[, 1] <- ifelse(exp_ci[, 1] < 0, 0, exp_ci[, 1])
    exp_ci[, 2] <- ifelse(exp_ci[, 2] > 1, 1, exp_ci[, 2])
    
    expect_equal(act_ci$CI, exp_ci, check.attribute = F)
    
  }
  
})


plot(r, act_ci$rec, type = "l")
lines(r, act_ci$CI[, 1], type = "l", col = "red")
lines(r, act_ci$CI[, 2], type = "l", col = "red")


plot(r, act_ci$prec, type = "l")
lines(r, act_ci$CI[, 1], type = "l", col = "red")
lines(r, act_ci$CI[, 2], type = "l", col = "red")
