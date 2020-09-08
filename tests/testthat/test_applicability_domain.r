context("applicability_domain")

test_that("Hotteling T2 outlier results matches results in paper", {
  
  skip_on_cran()
  data("aid364")
  train <- aid364[1:400, ]
  test.burd <- aid364[401:500, 3:26]
  
  cml <- ModelTrain(train, ids = T, xcol.lengths = c(24, 147),
                    des.names = c("BurdenNumbers", "Pharmacophores"), verbose = F)
  
  
  expect_doppelganger("applicability-domain", ApplicabilityDomain(traindata = cml$data[[1]],
                                                                  testdata = test.burd,
                                                                  desname = "Burden Numbers",
                                                                  pvalue = .01))
  outliers <- ApplicabilityDomain(traindata = cml$data[[1]],
                                  testdata = test.burd,
                                  desname = "Burden Numbers",
                                  pvalue = .01)
  
  act_outliers <- outliers$test.outliers
  exp_outliers <- c("404", "405", "409", "420", "421", "430", "437", "440", "448", 
                    "451", "460", "462", "465", "470", "471", "481", "487", "490", 
                    "492", "494", "495", "496", "497")
  
  expect_equal(act_outliers, exp_outliers)
})
