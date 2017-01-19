library(devtools)

packages <- c("tree", "rpart", "randomForest",
              "class", "e1071", "lars", "MASS",
              "elasticnet", "nnet", "pROC",
              "foreach")
packages <- packages[order(packages)]

for(pkg in packages){
  use_package(pkg)
}

for(pkg in packages){
  print(packageVersion(pkg))
}
