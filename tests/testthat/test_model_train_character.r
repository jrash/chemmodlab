context("ModelTrain with molecule input")

test_that("ModelTrain with molecule input matches the results in the paper", {
  
  skip_on_cran()
  # TODO make sure this works on your windows machine
  skip_on_travis()
  
  load("test_model_train_character.rdata")
  mols <- rcdk::parse.smiles(bpdata[, 1])
  bp <- bpdata[, 2]
  
  act_cml1 <- ModelTrain(descriptors = c("topological", "electronic"),
                    y = bp, mols = mols, verbose = F)
  
  act_cml2 <- ModelTrain(descriptors = c("fp.maccs", "fp.standard"),
                     y = bp, mols = mols, verbose = F)
  
  expect_equal(act_cml1, exp_cml1)
  expect_equal(act_cml2, exp_cml2)
})