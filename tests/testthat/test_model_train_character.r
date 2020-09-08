context("ModelTrain with molecule input")

test_that("ModelTrain with molecule input matches the results in the paper", {
  data(bpdata)
  mols <- rcdk::parse.smiles(bpdata[, 1])
  bp <- bpdata[, 2]
  
  cml1 <- ModelTrain(descriptors = c("topological", "electronic"),
                    y = bp, mols = mols)
  
  # Cleanup
  rm(bpdata)
})