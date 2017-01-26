
BackTree <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                     nperm, params) {
  #-----Recursive partitioning using "tree" with splitting criterion deviance and default settings. Specifically:
  # mincut=5, minimum leaf size
  # minsize=10, minimum parent size
  # mindev=.01, within-node deviance must be at least this times that of the root node for node to split
  # Possible modifications that have NOT been pursued here:
  # use cv.tree to perform k-fold CV for selecting tree size
  set.seed(current.seed)
  cat("Tree ---------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  #  temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  #  impdesc <- c()
  for (id in 1:nfolds) {
    if (classify == "Y") {
      # Fit classification tree
      work.meth <- tree::tree(as.factor(y) ~ ., data = work.data,
                        subset = (fold.id != id), method = "class")
      # class predictions are made for the data that was in the hold out fold Do not
      # know why the predicted classes can't be returned and then converted to numeric
      work.pred[fold.id == id] <-
        as.numeric(levels(as.factor(work.data[, 1])))[predict(work.meth,
                                                              work.data[fold.id == id, ],
                                                              type = "class")]
      # get the probability for the first class label (0)
      work.prob[fold.id == id] <- predict(work.meth, work.data[fold.id == id, ])[, 1]
    } else {
      # Fit regression tree
      work.meth <- tree::tree(y ~ ., data = work.data, subset = (fold.id != id))
      work.pred[fold.id == id] <- predict(work.meth,
                                          work.data[fold.id == id, ], type = "vector")
    }
  }
  work.model.acc <- BackAssess(work.data[,1], work.pred, work.impdesc, classify=classify)
  #returns probability of being active
  return(list(pred = work.pred, impdesc = work.impdesc,
              prob = 1 - work.prob, model.acc = work.model.acc))
}


#--------------------------------------------------------------------------------
BackRpart <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                      nperm, params) {
  #-----Recursive partitioning using "rpart" with splitting criterion "information" and
  # minbucket=5, minimum leaf size
  # minsplit=10, minimum parent size
  # maxcompete=0, don't get information on competitive splits
  # maxsurrogate=0, don't get information on surrogate splits
  # cp = .01,  must increase R^2 by .01 with each split
  # Possible modifications that have NOT been pursued here:
  # many ...
  set.seed(current.seed)
  cat("RPart --------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  for (id in 1:nfolds) {
    if (classify == "Y") {
      # setting maxsurrpgate to zero saves computational time
      work.meth <- rpart::rpart(as.factor(y) ~ ., data = work.data, subset = (fold.id != id),
                         method = "class", parms = list(split = "information"),
                         control = rpart::rpart.control(minsplit = 10,
                                                 minbucket = 5, maxcompete = 0,
                                                 maxsurrogate = 0, cp = params$RPart$cp))
      work.pred[fold.id == id] <-
        as.numeric(levels(as.factor(work.data[, 1])))[predict(work.meth,
                                                              work.data[fold.id == id, ], type = "class")]
      work.prob[fold.id == id] <- predict(work.meth, work.data[fold.id == id, ])[, 1]
    } else {
      work.meth <- rpart::rpart(y ~ ., data = work.data, subset = (fold.id != id),
                         method = "anova", control = rpart::rpart.control(minsplit = 10, minbucket = 5,
                                                                   maxcompete = 0, maxsurrogate = 0))
      work.pred[fold.id == id] <- predict(work.meth, work.data[fold.id == id, ], type = "vector")
    }
  }
  work.model.acc <- BackAssess( work.data[,1], work.pred, work.impdesc, classify=classify )
  return( list(pred=work.pred, impdesc=work.impdesc, prob=1-work.prob,  model.acc = work.model.acc) )
}

#--------------------------------------------------------------------------------
BackRf <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                   nperm, params) {
  #-----Random Forest using
  # ntree=100
  # mtry = sqrt(p) [for classification] and = p/3 [for regression]
  # nodesize = 5
  # importance = TRUE (to calculate important descriptors)
  # Possible modifications that have NOT been pursued here:
  # many ...
  set.seed(current.seed)
  cat("Forest -------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  for (id in 1:nfolds) {
    if (classify == "Y") {
      if (is.null(params$Forest$mtry)) {
        work.meth <-
          randomForest::randomForest(y=as.factor(work.data$y[fold.id!=id]),
                                  x=work.data[fold.id!=id, -1], ntree=100, nodesize=5)
      } else {
        work.meth <-
          randomForest::randomForest(y=as.factor(work.data$y[fold.id!=id]),
                                  x=work.data[fold.id!=id, -1], ntree=100, nodesize=5,
                                  mtry = params$Forest$mtry)
      }
      work.pred[fold.id==id] <-
        as.numeric(levels(as.factor(work.data$y)))[
          predict(work.meth, work.data[fold.id==id,-1])]
      work.prob[fold.id==id] <-
        predict(work.meth,
                              work.data[fold.id==id,-1], type="prob")[,1]
    } else {
      #TO DO: I dont think you need this if then anymore
      if (is.null(params$Forest$mtry)) {
        work.meth <- randomForest::randomForest(y=work.data$y[fold.id!=id],
                                  x=work.data[fold.id!=id, -1], ntree=100,
                                  nodesize=5)
      } else {
        work.meth <- randomForest::randomForest(y=work.data$y[fold.id!=id],
                                  x=work.data[fold.id!=id, -1], ntree=100,
                                  nodesize=5,
                                  mtry = params$Forest$mtry)
      }
      work.pred[fold.id==id] <- predict(work.meth,
                                                      work.data[fold.id==id,-1])
    }
  }
  work.model.acc <- BackAssess( work.data[,1], work.pred, work.impdesc, classify=classify )
  return( list(pred=work.pred, impdesc=work.impdesc,
               prob=1-work.prob,  model.acc = work.model.acc) )
}

#--------------------------------------------------------------------------------
BackSvm <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                    nperm, params) {
  #-----Support Vector Machine using
  # kernel = radial.basis
  # gamma = 1
  # cost = 1
  # epsilon = .01
  # Possible modifications that have NOT been pursued here:
  # many ...
  set.seed(current.seed)
  cat("SVM ----------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  # temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  for (id in 1:nfolds) {
    X <- work.data[fold.id != id, -1]
    # werent these already filtered out?
    varying.cols <- (apply(X, 2, var) > 0)
    if (classify == "Y") {
      work.meth <- e1071::svm(y = as.factor(work.data$y[fold.id != id]),
                              x = X[, varying.cols],
                       gamma = params$SVM$gamma,
                       cost = params$SVM$cost, probability = TRUE)
      # removing the response column and all non varying columns from descriptor matrix
      work.pred[fold.id == id] <-
        as.numeric(levels(as.factor(work.data$y)))[predict(work.meth,
                                                           work.data[fold.id == id,
                                                                     c(FALSE, varying.cols)])]
      temp.prob <- predict(work.meth, work.data[fold.id == id,
                                                       c(FALSE, varying.cols)],
                                  probability = TRUE)
      work.prob[fold.id == id] <- attr(temp.prob, "probabilities")[, 1]
    } else {
      work.meth <- e1071::svm(y = work.data$y[fold.id != id], x = X[, varying.cols],
                       gamma = 1, epsilon = params$SVM$epsilon)
      work.pred[fold.id == id] <- predict(work.meth, work.data[fold.id == id,
                                                               c(FALSE, varying.cols)])
    }
  }
  work.model.acc <- BackAssess( work.data[,1], work.pred, work.impdesc,
                                classify=classify )
  return( list(pred=work.pred, impdesc=work.impdesc,
               prob=1-work.prob,  model.acc = work.model.acc) )
}

#--------------------------------------------------------------------------------
BackNnet <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify,
                     current.seed, nperm, params) {
  #-----Neural Network using "nnet" with:
  # size=2
  # decay = 0
  # trace=FALSE (don't print out convergence info)
  # Possible modifications that have NOT been pursued here:
  # many
  # TODO, Why not for Regression?
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  # temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  if (classify == "Y") {
    cat("NNet ---------------------------------------", "\n")
    for (id in 1:nfolds) {
      X <- work.data[fold.id != id, -1]
      varying.cols <- (apply(X, 2, var) > 0)
      work.meth <- nnet::nnet(as.factor(y) ~ .,
                        data = work.data[fold.id != id, c(TRUE, varying.cols)],
                        size = params$NNet$size, decay = params$NNet$decay,
                        trace = FALSE)
      work.pred[fold.id == id] <- predict(work.meth,
                                          work.data[fold.id == id,
                                                    c(TRUE, varying.cols)],
                                          type = "class")
      work.prob[fold.id == id] <- 1 - predict(work.meth,
                                              work.data[fold.id == id,
                                                        c(TRUE, varying.cols)],
                                              type = "raw")
    }
    work.model.acc <- BackAssess(work.data[, 1], work.pred, work.impdesc,
                                 classify = classify)
  }
  return(list(pred = work.pred, impdesc = work.impdesc,
              prob = 1 - work.prob, model.acc = work.model.acc))
}


#--------------------------------------------------------------------------------
BackKnn <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                    nperm, params) {
  #-----K nearest neighbor using
  # k = 10
  # Possible modifications that have NOT been pursued here:
  # cv option ...
  # TODO why no knn regression?
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  # temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  if (classify == "Y") {
    cat("KNN ----------------------------------------", "\n")
    for (id in 1:nfolds) {
      work.meth <- class::knn(cl = as.factor(work.data$y[fold.id != id]),
                            train = work.data[fold.id != id, -1],
                       test = work.data[fold.id == id, -1], k = params$KNN$k,
                       prob = TRUE)
      work.pred[fold.id == id] <- as.numeric(levels(as.factor(work.data$y)))[work.meth]
      #If the predicted class is 1, use the probability, if it is 0, use 1- prob.
      #Gets the probability of class 1
      work.prob[fold.id==id] <- abs(work.pred[fold.id==id] +
                                       attr(work.meth, "prob") - 1)
    }

    work.model.acc <- BackAssess(work.data[, 1], work.pred, work.impdesc, classify = classify)
  }
  return(list(pred = work.pred, impdesc = work.impdesc,
              prob = work.prob, model.acc = work.model.acc))
}


#--------------------------------------------------------------------------------
BackLars <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                     nperm, params) {
  #-----Lars
  # max.steps = min(m,n-intercept)
  cat("LARs ---------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  # temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  for (id in 1:nfolds) {
    if (is.null(params$LARs)) {
      work.meth <- lars::lars(y = work.data[fold.id != id, 1],
                        x = as.matrix(work.data[fold.id != id, -1]), type = "lar")
    } else {
      work.meth <- lars::lars(y = work.data[fold.id != id, 1],
                        x = as.matrix(work.data[fold.id != id, -1]), type = "lar",
                        max.steps = params$LARs$steps)
    }

    temp.pred <- predict(work.meth, as.matrix(work.data[fold.id == id, -1]),
                         type = "fit")
    work.pred[fold.id == id] <- temp.pred$fit[, dim(temp.pred$fit)[2]]
  }
  work.model.acc <- BackAssess(work.data[, 1], work.pred, work.impdesc, classify = "N")
  return(list(pred = work.pred, impdesc = work.impdesc, model.acc = work.model.acc))
}

#--------------------------------------------------------------------------------
BackRidge <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                      nperm, params) {
  #-----Ridge Regression using
  # lambda = 0.1
  cat("Ridge --------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  for (id in 1:nfolds) {
    X <- work.data[fold.id != id, -1]
    varying.cols <- (apply(X, 2, var) > 0)
    work.meth <- MASS::lm.ridge(y ~ ., data = work.data[fold.id != id,
                                                  c(TRUE, varying.cols)],
                          lambda = params$Ridge$lambda)
    work.pred[fold.id == id] <- (as.matrix(work.data[fold.id == id,
                                                     c(FALSE, varying.cols)]) %*%
                                   (work.meth$coef/work.meth$scales)) +
      ((work.meth$ym - (work.meth$xm %*% (work.meth$coef/work.meth$scales)))[1, 1])
  }
  work.model.acc <- BackAssess(work.data[, 1], work.pred, work.impdesc, classify = "N")
  return(list(pred = work.pred, impdesc = work.impdesc, model.acc = work.model.acc))
}

#--------------------------------------------------------------------------------
BackEnet <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                     nperm, params) {
  # lambda = 1
  #-----enet
  cat("ENet ---------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  #  temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  for(id in 1:nfolds){
    X <- work.data[fold.id != id, -1]
    varying.cols <- (apply(X, 2, var) > 0)
    work.meth <- elasticnet::enet(y = work.data[fold.id != id, 1],
                                  x = as.matrix(X[, varying.cols]),
                                  lambda = params$ENet$lambda)
    temp.pred <- elasticnet::predict(work.meth,
                         as.matrix(work.data[fold.id == id,
                                             c(FALSE, varying.cols)]),
                         type = "fit")
    work.pred[fold.id == id] <- temp.pred$fit[, dim(temp.pred$fit)[2]]
  }
  work.model.acc <- BackAssess(work.data[, 1], work.pred, work.impdesc,
                               classify = "N")
  return(list(pred = work.pred, impdesc = work.impdesc,
              model.acc = work.model.acc))
}


#--------------------------------------------------------------------------------
BackPcr <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify,
                    current.seed, nperm, params) {
  #-----Principal components regression using home-grown code
  cat("PCR ----------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  #  temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  #  impdesc <- c()
  for (id in 1:nfolds) {
    X <- work.data[fold.id != id, -1]
    varying.cols <- (apply(X, 2, var) > 0)
    newX <- work.data[fold.id == id, -1]
    work.meth <- PcrZG(X = X[, varying.cols], Y = work.data[fold.id != id, 1],
                       newX = newX[, varying.cols])
    work.pred[fold.id == id] <- work.meth$Ypred
  }
  work.model.acc <- BackAssess(work.data[,1], work.pred, work.impdesc,
                                classify="N")
  return(list(pred=work.pred, impdesc=work.impdesc, model.acc = work.model.acc))
}

#--------------------------------------------------------------------------------
PcrZG <- function(X, Y, newX = NULL) {
  Xcs <- scale(X)
  Yc <- scale(Y, scale = FALSE)
  full.svd <<- La.svd(Xcs)
  nLV.ZG <- ZhuGhodsi((full.svd$d)^2)
  #  nLV.ZG <- sum( full.svd$d>0 )
  #  print("nLV.ZG"); print(nLV.ZG)
  U <- matrix(full.svd$u[, 1:nLV.ZG], nrow = nrow(full.svd$u), ncol = nLV.ZG)
  D <- as.vector(full.svd$d[1:nLV.ZG])
  Vt <- matrix(full.svd$vt[1:nLV.ZG, ], nrow = nLV.ZG, ncol = ncol(full.svd$vt))

  ScaleCen.B <- t(Vt) %*% diag(1/D, nrow = nLV.ZG) %*% t(U) %*% Yc
  original.B <- diag(1/attr(Xcs, "scaled:scale")) %*% ScaleCen.B

  df.error <- nrow(Xcs) - nLV.ZG - 1
  uut <- -U %*% t(U)
  for (i in 1:nrow(U)) uut[i, i] <- 1 + uut[i, i]
  mse <- t(Yc) %*% uut %*% Yc/df.error
  # mse <- t(Yc) %*% ( diag(nrow(Xcs)) - U%*%t(U) ) %*% Yc /df.error
  mse <- max(mse, 1e-16)
  seScaleCen.B <- sqrt(mse) * sqrt(diag(t(Vt) %*% diag(1/(D^2),
                                                       nrow = nLV.ZG) %*% Vt))

  singular.vals <- sqrt(apply((full.svd$u) %*% diag(full.svd$d), 2, var))
  imp.predictors <- 1 * (abs(ScaleCen.B) > qt(0.975, df.error) * seScaleCen.B)

  if (is.null(newX))
    Ypred <- attr(Yc, "scaled:center") +
    as.matrix(X - attr(Xcs, "scaled:center")) %*% original.B
  else
    Ypred <- attr(Yc, "scaled:center") +
    as.matrix(newX - attr(Xcs, "scaled:center")) %*% original.B

  return(list(Ypred = Ypred, imp.predictors = imp.predictors))
}


#--------------------------------------------------------------------------------
ZhuGhodsi <- function(x) {
  #-----Zhu & Ghodsi method for determining the number of latent variables given
  # some measure of relative importance in x, where x1>=x2>=...
  # Computational Statistics and Data Analysis 2006?
  n <- length(x)
  nLV <- 1:n
  pVar.q <- rep(NA, n)
  pVar.q[1] <- var(x)  #q=0
  pVar.q[2] <- var(x[2:n])  #q=1
  pVar.q[n] <- var(x[1:(n - 1)])  #q=n-1
  for (q in 2:(n - 2)) {
    i <- q + 1
    pVar.q[i] <- ((q - 1) * var(x[1:q]) +
                    (n - q - 1) * var(x[(q + 1):n]))/(n - 2)
  }
  #  print(pVar.q)
  #  return( (1:n)[(order(pVar.q))[1]] )
  return(-1 + (order(pVar.q))[1])
}


#--------------------------------------------------------------------------------
BackPlsR <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                     nperm, params) {
  #-----Partial least squares using home-grown code based on 'kernelpls'
  cat("PLS ----------------------------------------", "\n")
  work.pred <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  #  temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  #  temp.impdesc <- matrix(0,n.pred,nfolds)
  for (id in 1:nfolds) {
    X <- work.data[fold.id != id, -1]
    varying.cols <- (apply(X, 2, var) > 0)
    newX <- work.data[fold.id == id, -1]
    work.meth <- KernelPlsNew(X = X[, varying.cols],
                              Y = work.data[fold.id != id, 1],
                              ncomp = min(n.obs, n.pred, 100),
                              newX = newX[, varying.cols])
    #      work.pred[fold.id==id] <- work.meth$Ypred[,1,min(n.obs,n.pred,100)]
    #      work.pred[fold.id==id] <- work.meth$Ypred
    nLV.ZG <- ZhuGhodsi(work.meth$Yvar)
    #      plot( work.meth$Yvar / work.meth$Ytotvar ); abline(v=nLV.ZG)
    work.pred[fold.id==id] <- work.meth$Ypred[, 1, nLV.ZG]
  }
  work.model.acc <- BackAssess(work.data[,1], work.pred,
                               work.impdesc, classify="N")
  #  plot( work.data[,1], work.pred )
  return(list(pred=work.pred, impdesc=work.impdesc, model.acc = work.model.acc))
}

#--------------------------------------------------------------------------------
KernelPlsNew <- function(X, Y, ncomp, newX, stripped = FALSE, ...) {
  # Based on Algorithm 1 of Dayal and MacGregor 1997 J of Chemometrics 11, 73-85

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  newX <- as.matrix(newX)
  if (!stripped) {
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
  }
  dimnames(X) <- dimnames(Y) <- NULL
  nobj <- dim(X)[1]
  nvar <- dim(X)[2]
  npred <- dim(Y)[2]
  X <- scale(X)
  Xmeans <- attr(X, "scaled:center")
  Xscales <- attr(X, "scaled:scale")
  Y <- scale(Y)
  Ymeans <- attr(Y, "scaled:center")
  Yscales <- attr(Y, "scaled:scale")
  PP <- matrix(0, ncol = ncomp, nrow = nvar)
  QQ <- matrix(0, ncol = ncomp, nrow = npred)
  TT <- matrix(0, ncol = ncomp, nrow = nobj)
  RR <- matrix(0, ncol = ncomp, nrow = nvar)
  B <- array(0, c(dim(X)[2], dim(Y)[2], ncomp))
  Ypred <- array(0, c(dim(newX)[1], npred, ncomp))
  if (!stripped)
    UU <- matrix(0, ncol = ncomp, nrow = nobj)
  XtY <- crossprod(X, Y)
  XtYtotvar <- sum((XtY)^2)
  for (a in 1:ncomp) {
    if (npred == 1)
      ww <- XtY else {
        ww <- La.svd(XtY)$u[, 1]
        #            if (npred > nvar)
        #                qq <- La.svd(XtY)$vt[1, ]
        #            else qq <- La.svd(XtY)$u[1, ]
        #            ww <- XtY %*% qq
      }
    rr <- ww
    if (a > 1)
      for (j in 1:(a - 1)) rr <- rr - (sum(PP[, j] * ww)) * RR[, j]
    tt <- X %*% rr
    ttnorm <- (sum(tt * tt))
    pp <- crossprod(X, tt)/ttnorm
    qq <- crossprod(XtY, rr)/ttnorm
    XtY <- XtY - (pp %*% t(qq)) * ttnorm
    uu <- Y %*% matrix(qq, ncol = 1)  # I'm not sure about this line!
    if (a > 1)  # Im not sure about this line!
      uu <- uu - TT %*% crossprod(TT, uu)  # Im not sure about this  ine!
    TT[, a] <- tt
    PP[, a] <- pp
    QQ[, a] <- qq
    RR[, a] <- rr
    if (!stripped)
      UU[, a] <- uu
    B[, , a] <- RR[, 1:a, drop = FALSE] %*% t(QQ[, 1:a, drop = FALSE])
    temp.Ypred.a <- as.matrix(scale(newX, center = Xmeans, scale = Xscales) %*%
                                B[, , a])
    temp.Ypred.a <- sweep(temp.Ypred.a, 2, Yscales, "*")
    Ypred[, , a] <- sweep(temp.Ypred.a, 2, Ymeans, "+")
  }
  if (stripped) {
    list(coefficients = B, Ypred = Ypred)
  } else {
    objnames <- dnX[[1]]
    if (is.null(objnames))
      objnames <- dnY[[1]]
    xvarnames <- dnX[[2]]
    yvarnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(TT) <- dimnames(UU) <- list(objnames, compnames)
    dimnames(RR) <- dimnames(PP) <- list(xvarnames, compnames)
    dimnames(QQ) <- list(yvarnames, compnames)
    dimnames(B) <- list(xvarnames, yvarnames, nCompnames)
    dimnames(Ypred) <- list(NULL, yvarnames, nCompnames)
    class(TT) <- class(UU) <- "scores"
    class(PP) <- class(QQ) <- "loadings"
    colSumsTT <- colSums(TT^2)
    list(coefficients = B, Ypred = Ypred, Xmeans = Xmeans, Ymeans = Ymeans,
         Xscores = TT, Xloadings = PP, Yscores = UU, Yloadings = QQ, projection = RR,
         Xvar = colSums(PP^2)*colSumsTT, Xtotvar = sum(X^2),
         Yvar = colSums(QQ^2)*colSumsTT, Ytotvar = sum(Y^2) )
  }
}


#--------------------------------------------------------------------------------
BackPlsLdaNew <- function(work.data, n.obs, n.pred, nfolds, fold.id, classify, current.seed,
                          nperm, params) {
  #-----Partial least squares using home-grown code based on 'kernelpls'
  work.pred <- rep(NA, n.obs)
  work.prob <- rep(NA, n.obs)
  work.impdesc <- c()
  set.seed(current.seed)
  #  temp.impdesc <- matrix( rep(NA,n.pred*nperm*nfolds), ncol=n.pred )
  #  temp.impdesc <- matrix(0,n.pred,nfolds)
  if (classify == "Y") {
    cat("PLSLDA -------------------------------------", "\n")
    for (id in 1:nfolds) {
      X <- work.data[fold.id != id, -1]
      varying.cols <- (apply(X, 2, var) > 0)
      newX <- work.data[fold.id == id, -1]
      X.train <- X[, varying.cols]
      Y.train <- work.data[fold.id != id, 1]
      newX.test <- newX[, varying.cols]
      work.meth <- SimPlsLdaNew(X = X.train, Y = Y.train,
                                ncomp = min(n.obs, n.pred, 100),
                                newX = newX.test)
      nLV.ZG <- max(ZhuGhodsi(work.meth$Yvar), 1)
      work.prob[fold.id==id] <- work.meth$Yprob[, 2, nLV.ZG]
      work.pred[fold.id==id] <- work.meth$Ypred[, nLV.ZG]
    }
    work.model.acc <- BackAssess(work.data[, 1], work.pred, work.impdesc,
                                 classify="Y")
  }
  return(list(pred=work.pred, impdesc=work.impdesc, prob=work.prob,
              model.acc = work.model.acc))
}


#--------------------------------------------------------------------------------
SimPlsLdaNew <- function(X, Y, ncomp, newX, stripped = FALSE, priors = NULL, ...) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Y.orig <- Y
  newX <- as.matrix(newX)

  # Get information on classes. This code assumes interest is only in Y[,1] and
  # this single column may indicate many different classes
  class.labels <- unique(Y[, 1])
  num.classes <- length(class.labels)
  # id.in.class <- matrix( 0, nrow(Y), num.classes )
  id.in.class <- array(NA, c(nrow(Y), num.classes))
  for (j in 1:num.classes) {
    id.in.class[, j] <- Y[, 1] == class.labels[j]
  }
  num.in.classes <- colSums(id.in.class)

  # Initiate LDA output by determining priors
  if (is.null(priors)) {
    priors <- num.in.classes/sum(num.in.classes)
    logpriors <- log(priors)
  } else {
    logpriors <- log(priors)
  }

  if (!stripped) {
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
  }
  dimnames(X) <- dimnames(Y) <- NULL
  nobj <- dim(X)[1]
  nvar <- dim(X)[2]
  npred <- dim(Y)[2]
  ntest <- dim(newX)[1]
  X <- scale(X)
  Xmeans <- attr(X, "scaled:center")
  Xscales <- attr(X, "scaled:scale")
  ScaleCenter.newX <- as.matrix(scale(newX, center = Xmeans, scale = Xscales))
  Y <- scale(Y)
  Ymeans <- attr(Y, "scaled:center")
  Yscales <- attr(Y, "scaled:scale")
  PP <- matrix(0, ncol = ncomp, nrow = nvar)
  QQ <- matrix(0, ncol = ncomp, nrow = npred)
  TT <- matrix(0, ncol = ncomp, nrow = nobj)
  RR <- matrix(0, ncol = ncomp, nrow = nvar)
  Yprob <- array(0, c(ntest, num.classes, ncomp))
  Ypred <- matrix(0, ntest, ncomp)
  if (!stripped)
    UU <- matrix(0, ncol = ncomp, nrow = nobj)
  XtY <- crossprod(X, Y)
  XtYtotvar <- sum((XtY)^2)

  simpls.scores <- pls::simpls.fit(X, Y, ncomp)
  TT <- simpls.scores$scores
  RR <- simpls.scores$projection
  newTT <- ScaleCenter.newX %*% as.matrix(RR)
  simpls.lda <- MASS::lda(x = TT, group = Y.orig)
  for (a in 1:ncomp) {
    lda.a <- predict(object = simpls.lda, newdata = newTT, dimen = a)
    Yprob[, , a] <- lda.a$posterior
    Ypred[, a] <- lda.a$class
  }
  if (stripped) {
    list(Yprob = Yprob, Ypred = Ypred)
  } else {
    objnames <- dnX[[1]]
    if (is.null(objnames))
      objnames <- dnY[[1]]
    xvarnames <- dnX[[2]]
    yvarnames <- dnY[[2]]
    compnames <- paste("Comp", 1:ncomp)
    nCompnames <- paste(1:ncomp, "comps")
    dimnames(TT) <- dimnames(UU) <- list(objnames, compnames)
    dimnames(RR) <- dimnames(PP) <- list(xvarnames, compnames)
    dimnames(QQ) <- list(yvarnames, compnames)
    dimnames(Ypred) <- list(NULL, nCompnames)
    class(TT) <- class(UU) <- "scores"
    class(PP) <- class(QQ) <- "loadings"
    colSumsTT <- colSums(TT^2)
    list(class.labels = class.labels, Yprob = Yprob, Ypred = Ypred,
         Xmeans = Xmeans, Ymeans = Ymeans, Xscores = TT, Xloadings = PP,
         Yscores = UU,
         Yloadings = QQ, projection = RR, Xvar = colSums(PP^2)*colSumsTT,
         Xtotvar = sum(X^2),
         Yvar = colSums(QQ^2)*colSumsTT, Ytotvar = sum(Y^2))
  }
}
