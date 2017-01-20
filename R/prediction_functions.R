# Disregarding prediction at the moment

BackPredict <- function(yfilein, ycol, xfilein, xcols, filepred, fileimpdesc, fileprob,
                        infofile = NA, idcol = NA, modelyfile, modelxfile) {
  #-----Background predicting program
  # yfilein=input response data: response in column ycol xfilein=input descriptors
  # data: descriptors in columns xcols filepred=output file of predictions:
  # observed response in column 1, one additional column of predictions for each
  # method modelyfile=model response data modelxfile=model descriptors data
  # idcol=column which contains compound labels

  #-----Read in methods
  if (missing(infofile))
    meths <- c("NNet", "PCR", "ENet", "PLS", "Ridge", "LARs", "PLSLDA", "RPart",
               "Tree", "SVM", "KNN", "Forest", "Forest70", "TreeEns", "RPartEns", "KNNEns")
  else meths <- read.table(file = infofile, sep = "|", skip = 1, nrows = 1, colClasses = "character")

  #-----Read in model data
  model.list <- read.in.data(yfilein = modelyfile, xfilein = modelxfile, idcol = idcol)
  model.data <- model.list[[1]]

  #-----Read in prediction data, and discard constant columns from model data
  work.data <- read.in.data(yfilein, ycol, xfilein, xcols, idcol = idcol, type = "PREDICT")[[1]]
  work.data <- subset(work.data, select = (!model.list[[2]]))

  n.obs <- nrow(work.data)
  n.pred <- ncol(work.data) - 1



  #-----Determine if we have binary or continuous data
  if (sum(!(model.data$y %in% c(1, 0))))
    classify <- "N" else classify <- "Y"

  #-----Start outputing to summary file
  cat("Begining Prdeictions: ", date(), "\n")
  cat("Model Descriptors: ", modelxfile, "(", (ncol(model.data) - 1), ")\n")
  cat("Model Responses: ", modelyfile, "(", nrow(model.data), ")\n")
  cat("Prediction Descriptors: ", xfilein, "(", n.pred, ")\n")
  if ((!missing(yfilein)) && (yfilein != ""))
    cat("Prediction Responses: ", yfilein, "(", n.obs, ")\n\n\n")
  all.preds <- data.frame(Observed = work.data$y)
  if (classify == "Y")
    all.probs <- data.frame(Observed = work.data$y) else all.probs <- NA

  IDS <- rownames(work.data)

  #-----'tree' method
  if (sum(meths %in% "Tree") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredTree(model.data, work.data,
                                                        classify), error = function(e) {
                                                          warning(paste("WARNING...Tree not run:", e$message))
                                                          work.results <- list()
                                                        }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, Tree = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, Tree = work.results$prob)
    }
  }

  #-----'rpart' method
  if (sum(meths %in% "RPart") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredRpartCml(model.data, work.data,
                                                            classify), error = function(e) {
                                                              warning(paste("WARNING...RPart not run:", e$message))
                                                              work.results <- list()
                                                            }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, RPart = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, RPart = work.results$prob)
    }
  }

  #-----'randomforest' method
  if (sum(meths %in% "Forest") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredRf(model.data, work.data,
                                                      classify), error = function(e) {
                                                        warning(paste("WARNING...RF not run:", e$message))
                                                        work.results <- list()
                                                      }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, RF = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, RF = work.results$prob)
    }
  }

  #-----'svm' method
  if (sum(meths %in% "SVM") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredSvm(model.data, work.data,
                                                       classify), error = function(e) {
                                                         warning(paste("WARNING...SVM not run:", e$message))
                                                         work.results <- list()
                                                       }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, SVM = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, SVM = work.results$prob)
    }
  }

  #-----'nnet' method
  if ((sum(meths %in% "NNet") == 1) && (classify == "Y")) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredNnet(model.data, work.data,
                                                        classify), error = function(e) {
                                                          warning(paste("WARNING...NNet not run:", e$message))
                                                          work.results <- list()
                                                        }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, NNet = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, NNet = work.results$prob)
    }
  }

  #-----'knn' method
  if ((sum(meths %in% "KNN") == 1) && (classify == "Y")) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredKnn(model.data, work.data,
                                                       classify), error = function(e) {
                                                         warning(paste("WARNING...KNN not run:", e$message))
                                                         work.results <- list()
                                                       }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, KNN = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, KNN = work.results$prob)
    }
  }

  #-----'pls.lda' method
  if ((sum(meths %in% "PLSLDA") == 1) && (classify == "Y")) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredPlsLdaNew(model.data, work.data,
                                                             classify), error = function(e) {
                                                               warning(paste("WARNING...PLSLDA not run:", e$message))
                                                               work.results <- list()
                                                             }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, PLSLDA = work.results$pred)
      if (classify == "Y")
        all.probs <- data.frame(all.probs, PLSLDA = work.results$prob)
    }
  }

  #-----'lars' method
  if (sum(meths %in% "LARs") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredLars(model.data, work.data,
                                                        classify), error = function(e) {
                                                          warning(paste("WARNING...LAR not run:", e$message))
                                                          work.results <- list()
                                                        }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, LAR = work.results$pred)
    }
  }

  #-----'lm.ridge' method
  if (sum(meths %in% "Ridge") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredRidge(model.data, work.data,
                                                         classify), error = function(e) {
                                                           warning(paste("WARNING...Ridge not run:", e$message))
                                                           work.results <- list()
                                                         }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, Ridge = work.results$pred)
    }
  }

  #-----'enet' method
  if (sum(meths %in% "ENet") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredEnet(model.data, work.data,
                                                        classify), error = function(e) {
                                                          warning(paste("WARNING...ENet not run:", e$message))
                                                          work.results <- list()
                                                        }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, ENet = work.results$pred)
    }
  }

  #-----'PcrZG' method
  if (sum(meths %in% "PCR") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredPcr(model.data, work.data,
                                                       classify), error = function(e) {
                                                         warning(paste("WARNING...PCR not run:", e$message))
                                                         work.results <- list()
                                                       }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, PCR = work.results$pred)
    }
  }

  #-----'pls.R' method
  if (sum(meths %in% "PLS") == 1) {
    work.results <- list()
    pt <- proc.time()
    st <- system.time(tryCatch(work.results <- PredPlsR(model.data, work.data,
                                                        classify), error = function(e) {
                                                          warning(paste("WARNING...PLS not run:", e$message))
                                                          work.results <- list()
                                                        }))
    PrintTime(pt, st)
    if (length(work.results) > 0) {
      all.preds <- data.frame(all.preds, PLS = work.results$pred)
    }
  }

  cat("Ending Predictions: ", date(), "\n\n\n")
  #-----Write 'predictions', 'probabilities' and 'important descriptors' files
  desc <- xfilein
  while ((pos <- regexpr("/", desc, fixed = TRUE)) > 0) desc <- substr(desc, (pos +
                                                                                1), nchar(desc))
  desc <- substr(desc, 1, (regexpr(".", desc, fixed = TRUE) - 1))
  head <- paste(",,", desc, sep = "")
  if (ncol(all.preds) > 2)
    for (i in 2:(ncol(all.preds) - 1)) head <- paste(head, ",", sep = "")
  write.table(head, file = filepred, quote = FALSE, col.names = FALSE, row.names = FALSE)
  all.preds <- cbind(IDS, all.preds)
  write.table(all.preds, file = filepred, quote = FALSE, sep = ",", row.names = FALSE,
              append = TRUE)
  if (classify == "Y") {
    head <- paste(",,", desc, sep = "")
    if (ncol(all.probs) > 2)
      for (i in 2:(ncol(all.probs) - 1)) head <- paste(head, ",", sep = "")
      write.table(head, file = fileprob, quote = FALSE, col.names = FALSE, row.names = FALSE)
      all.probs <- cbind(IDS, all.probs)
      write.table(all.probs, file = fileprob, quote = FALSE, sep = ",", row.names = FALSE,
                  append = TRUE)
  }

  return(list(all.preds = all.preds, all.probs = all.probs, classify = classify,
              responses = work.data$y))
}


PredTree <- function(model.data, work.data, classify) {
  #-----Recursive partitioning using 'tree' with splitting criterion deviance and default settings. Specifically:
  # mincut=5, minimum leaf size minsize=10, minimum parent size mindev=.01,
  # within-node deviance must be at least this times that of the root node for node
  # to split
  cat("Tree ---------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- tree(as.factor(y) ~ ., data = model.data, method = "class")
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
                                                                     work.data, type = "class")]
    work.prob <- predict(work.meth, work.data)[, 1]
  } else {
    work.meth <- tree(y ~ ., data = work.data)
    work.pred <- predict(work.meth, work.data, type = "vector")
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}

#--- RPart has an internal function named 'PredRpart', our function is PredRpart.chm
PredRpartCml <- function(model.data, work.data, classify) {
  #-----Recursive partitioning using 'rpart' with splitting criterion 'information' and
  # minbucket=5, minimum leaf size minsplit=10, minimum parent size maxcompete=0,
  # don't get information on competitive splits maxsurrogate=0, don't get
  # information on surrogate splits Possible modifications that have NOT been
  # pursued here: many ...
  cat("RPart --------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- rpart(as.factor(y) ~ ., data = model.data, method = "class",
                       parms = list(split = "information"), minsplit = 10, minbucket = 5, maxcompete = 0,
                       maxsurrogate = 0)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
                                                                     work.data, type = "class")]
    work.prob <- predict(work.meth, work.data)[, 1]
  } else {
    work.meth <- rpart(y ~ ., data = model.data, method = "anova", parms = list(split = "information"),
                       minsplit = 10, minbucket = 5, maxcompete = 0, maxsurrogate = 0)
    work.pred <- predict(work.meth, work.data, type = "vector")
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredRf <- function(model.data, work.data, classify) {
  #-----Random Forest using
  # ntree=100 mtry = sqrt(p) [for classification] and = p/3 [for regression]
  # nodesize = 5 importance = TRUE (to calculate important descriptors) Possible
  # modifications that have NOT been pursued here: many ...
  cat("Forest -------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- randomForest(y = as.factor(model.data$y), x = model.data[, -1],
                              ntree = 100, nodesize = 5)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
                                                                     work.data[, -1])]
    work.prob <- predict(work.meth, work.data[, -1], type = "prob")[, 1]
  } else {
    work.meth <- randomForest(y = model.data$y, x = model.data[, -1], ntree = 100,
                              nodesize = 5)
    work.pred <- predict(work.meth, work.data[, -1])
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredSvm <- function(model.data, work.data, classify) {
  #-----Support Vector Machine using
  # kernel = radial.basis gamma = 1 Possible modifications that have NOT been
  # pursued here: many ...
  cat("SVM ----------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- svm(y = as.factor(model.data$y), x = model.data[, -1], gamma = 1,
                     probability = TRUE)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict(work.meth,
                                                                     work.data[, -1])]
    temp.prob <- predict(work.meth, work.data[, -1], probability = TRUE)
    work.prob <- attr(temp.prob, "probabilities")[, 1]
  } else {
    work.meth <- svm(y = model.data$y, x = model.data[, -1], gamma = 1)
    work.pred <- predict(work.meth, work.data[, -1])
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredNnet <- function(model.data, work.data, classify) {
  #-----Neural Network using 'nnet' with:
  # size=2 trace=FALSE (don't print out convergence info) Possible modifications
  # that have NOT been pursued here: many
  cat("NNet ---------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.pred <- rep(NA, n.work)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- nnet(as.factor(y) ~ ., data = model.data, size = 2, trace = FALSE)
    work.pred <- predict(work.meth, work.data, type = "class")
    work.prob <- 1 - predict(work.meth, work.data, type = "raw")
  }
  return(list(pred = work.pred, prob = 1 - work.prob))
}


PredKnn <- function(model.data, work.data, classify) {
  #-----K nearest neighbor using
  # k = 3 Possible modifications that have NOT been pursued here: cv option ...
  cat("KNN ----------------------------------------", "\n")
  n.work <- nrow(work.data)
  work.pred <- rep(NA, n.work)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.meth <- knn(cl = as.factor(model.data$y), train = model.data[, -1],
                     test = work.data[, -1], k = 10, prob = TRUE)
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[work.meth]
    work.prob <- abs(work.pred + attr(work.meth, "prob") - 1)
  }
  return(list(pred = work.pred, prob = work.prob))
}


PredKnnflex <- function(model.data, work.data, classify) {
  #-----K nearest neighbor using
  # k = 3
  cat("KNN ----------------------------------------", "\n")
  work.meth <- knn.dist(rbind(model.data[, -1], work.data[, -1]))
  n.model <- nrow(model.data)
  n.work <- nrow(work.data)
  work.prob <- rep(NA, n.work)
  if (classify == "Y") {
    work.pred <- as.numeric(knn.predict(train = 1:n.model, test = (n.model +
                                                                     1):(n.work + n.model), y = work.data$y, k = 10, dist.matrix = work.meth,
                                        agg.meth = "majority"))
    work.prob <- knn.predict(train = 1:n.model, test = (n.model + 1):(n.work +
                                                                        n.model), y = work.data$y, k = 10, dist.matrix = work.meth, agg.meth = "mean")
  } else {
    work.pred <- knn.predict(train = 1:n.model, test = (n.model + 1):(n.work +
                                                                        n.model), y = work.data$y, k = 10, dist.matrix = work.meth, agg.meth = "mean")
  }
  return(list(pred = work.pred, prob = work.prob))
}


PredLars <- function(model.data, work.data, classify) {
  #-----Lars
  cat("LARs ---------------------------------------", "\n")
  work.meth <- lars(y = model.data[, 1], x = as.matrix(model.data[, -1]), type = "lar")
  temp.pred <- predict(work.meth, as.matrix(work.data[, -1]), type = "fit")
  work.pred <- temp.pred$fit[, dim(temp.pred$fit)[2]]
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredRidge <- function(model.data, work.data, classify) {
  #-----Ridge Regression using
  # lambda = 0.1
  cat("Ridge --------------------------------------", "\n")
  work.meth <- lm.ridge(y ~ ., data = model.data, lambda = 0.1)
  work.pred <- (as.matrix(work.data[, -1]) %*% (work.meth$coef/work.meth$scales)) +
    ((work.meth$ym - (work.meth$xm %*% (work.meth$coef/work.meth$scales)))[1,
                                                                           1])
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredEnet <- function(model.data, work.data, classify) {
  #-----enet
  cat("ENet ---------------------------------------", "\n")
  work.meth <- enet(y = model.data$y, x = as.matrix(model.data[, -1]), lambda = 1)
  temp.pred <- predict(work.meth, as.matrix(work.data[, -1]), type = "fit")
  work.pred <- temp.pred$fit[, dim(temp.pred$fit)[2]]
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredPcr <- function(model.data, work.data, classify) {
  #-----Principal components regression using home-grown code
  cat("PCR ----------------------------------------", "\n")
  newX <- model.data[, -1]
  work.meth <- PcrZG(X = model.data[, -1], Y = model.data$y, newX = work.data[,
                                                                              -1])
  work.pred <- work.meth$Ypred
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredPlsR <- function(model.data, work.data, classify) {
  #-----Partial least squares using home-grown code based on 'kernelpls'
  cat("PLS ----------------------------------------", "\n")
  work.meth <- KernelPlsNew(X = model.data[, -1], Y = model.data$y, ncomp = min(nrow(model.data),
                                                                                (ncol(model.data) - 1), 100), newX = work.data[, -1])
  nLV.ZG <- ZhuGhodsi(work.meth$Yvar)
  work.pred <- work.meth$Ypred[, 1, nLV.ZG]
  return(list(pred = work.pred))
}


#--------------------------------------------------------------------------------
PredPlsLdaNew <- function(model.data, work.data, classify) {
  #-----Partial least squares using home-grown code based on 'kernelpls'
  cat("PLSLDA -------------------------------------", "\n")
  work.meth <- SimPlsLdaNew(X = model.data[, -1], Y = model.data$y, ncomp = min(nrow(model.data),
                                                                                (ncol(model.data) - 1), 100), newX = work.data[, -1])
  nLV.ZG <- max(ZhuGhodsi(work.meth$Yvar), 1)
  work.prob <- work.meth$Yprob[, 2, nLV.ZG]
  work.pred <- work.meth$Ypred[, nLV.ZG]
  return(list(pred = work.pred, prob = work.prob))
}
