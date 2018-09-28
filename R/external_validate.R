#' @import methods
#' @import stats
#' @import utils
#' 
#' @export
ExtValidate <- function(...) UseMethod("ExtValidate")

#' @describeIn ModelTrain Default S3 method
#' @export
ExtValidate.default <- function(x, cml,
                               models = c("NNet", "PLS", "LAR", "Lasso",
                                          "PLSLDA", "Tree", "SVM", "KNN", "RF"),
                               ...) {
  
  s3method <- "default"
  # checking parameters specified correctly
  # TODO check if each column is numeric in each matrix?
  # or should we keep just converting to numeric matrix?
  if (!is(x, "list")) {
    stop("'x' must be a list of numeric matrices")
  } else {
    for (i in seq_along(x)) {
      if (!is(x[[i]], "matrix"))
        stop("'x' must be a list of numeric matrices")
    }
  }
  
  BackExtValidate(x = x, cml = cml, models = models,
                 s3method = s3method)
}

#' @describeIn ModelTrain S3 method for class 'data.frame'
#' @export
ExtValidate.data.frame <- function(d,
                                  ids = FALSE,
                                  xcol.lengths = ifelse(ids,
                                                        length(d) - 2,
                                                        length(d) - 1),
                                  xcols = NA,
                                  models = c("NNet", "PLS", "LAR", "Lasso",
                                             "PLSLDA", "Tree", "SVM", "KNN", "RF"),
                                  ...) {
  
  #TODO dont force user to supply xcols
  s3method <- "data.frame"
  # checking parameters specified correctly
  if (!is(d, "data.frame"))
    stop("'d' must be a data frame")
  if (!is(xcol.lengths, "vector"))
    stop("'xcol.lengths' must be a vector of integers") 
  else {
    for (i in 1:length(xcol.lengths)) {
      if (!(xcol.lengths[[i]]%%1 == 0) || !is.numeric(xcol.lengths[[i]]))
        stop("'xcol.lengths' must be a list of integers")
    }
  }
  if (!is.na(xcols[[1]])) {
    if (!is(xcols, "list")) {
      stop("'xcols' must be a list of integer vectors") 
    } else {
      for (i in 1:length(xcols)) {
        if (!(all.equal(xcols[[i]], as.integer(xcols[[i]]))))
          stop("'xcols' must be a list of integer vectors")
      }
    }
  }
  if (!is(ids, "logical")) {
    stop("'ids' should be a logical")
  }
  if (ids == F) {
    if (sum(c(1, xcol.lengths)) > ncol(d))
      stop("there number of columns given is larger than number of columns in 'd'")
  } else {
    if (sum(c(2, xcol.lengths) > ncol(d)))
      stop("there number of columns given is larger than number of columns in 'd'")
  }
  
  idcol <- ifelse(ids, 1, NA)
  
  if (is.na(xcols)) {
    # use the descriptor column numbers to find the corresponding columns in the 
    # data frame
    xcols <- list()
    xcols[[1]] <- 1:xcol.lengths[1]
    if (length(xcol.lengths) > 1){
      for(i in 2:length(xcol.lengths)){
        l1 <- xcols[[i-1]][xcol.lengths[i-1]]
        l2 <- xcol.lengths[i]
        xcols[[i]] <- (l1+1):(l1+l2)
      }
    }
    # shift up indices by 1 if there is a id column
    for(i in 1:length(xcols)) {
      if(is.na(idcol)) {
        xcols[[i]] <- xcols[[i]]
      } else {
        xcols[[i]] <- xcols[[i]] + 1
      }
    }
  }
  BackExtValidate(d = d, ids = ids, cml = cml, 
                  xcol.lengths = xcol.lengths, xcols = xcols,
                  models = models, s3method = s3method, idcol = idcol)
}

BackExtValidate <- function(x = NA, d = NA, ids = NA, cml,
                           xcol.lengths = NA, xcols = NA,
                           idcol = NA, models, s3method) {
  
  if (s3method == "data.frame") {
    n.des <- length(xcols)
  } else {
    n.des <- length(x)
  }
  # checking parameters specified correctly
  if (!all(models %in% c("NNet", "PCR", "ENet", "PLS", "Ridge", "LAR", "PLSLDA",
                         "RPart", "Tree", "SVM", "KNN", "RF", "Lasso"))) {
    stop("'models' should be a character vector containing models existing in chemmodlab")
  }
  
  if (s3method == "data.frame") {
    n.obs <- nrow(d)
  } else {
    n.obs <- nrow(x[[1]])
  }
  
  #-----Determine if we have binary or continuous response
  # if sum == 0 then all of the response is either 0 or 1.  
  if (!exists("classify")) {
    if (sum(!(cml$responses %in% c(1, 0))))
      classify <- F 
    else 
      classify <- T
  }
  
  des.preds.ls <- list()
  des.probs.ls <- list()
    
  for (des.idx in 1:n.des) {
    # take response and current descriptor set columns
    if (s3method == "data.frame") {
      work.data <- ReadInData(d, ycol = NA, xcols[[des.idx]], idcol)[[1]]
    } else {
      work.data <- data.frame(x[[des.idx]])
    }
    n.pred <- ncol(work.data)
    
    #-----Make model parameter list
    params <- cml$params
    
    all.preds <- matrix(ncol = length(models), nrow = n.obs)
    colnames(all.preds) <- models 
    all.preds <- as.data.frame(all.preds)
    if (classify) {
      all.probs <- matrix(ncol = length(models), nrow = n.obs)
      colnames(all.probs) <- models
      all.probs <- as.data.frame(all.probs)
    }
    else all.probs <- NA
    
    IDS <- rownames(work.data)
    
    #### This section is calling the machine learning functions and saving the results.
    #### If results are succesfully returned, the predictions are appended to the
    #### prediction data frame
    
    model.x <- cml$data[[des.idx]]
    model.y <- cml$responses
    #-----'tree' method
    if (sum(models %in% "Tree") == 1) {
      work.results <- list()
      st <- tryCatch(work.results <- PredTree(model.y, model.x, work.data,
                                              classify, params),
                                 error = function(e) {
                                   warning(paste("WARNING...Tree not run:", e$message))
                                   work.results <- list()
                                 })
      if (length(work.results) > 0) {
        all.preds$Tree <- work.results$pred
        if (classify)
          all.probs$Tree <- work.results$prob
      }
    }
    if (sum(models %in% "RPart") == 1) {
      work.results <- list()
      st <- tryCatch(work.results <- PredRpart(model.y, model.x, work.data,
                                              classify, params),
                     error = function(e) {
                       warning(paste("WARNING...RPart not run:", e$message))
                       work.results <- list()
                     })
      if (length(work.results) > 0) {
        all.preds$RPart <- work.results$pred
        if (classify)
          all.probs$RPart <- work.results$prob
      }
    }
    
    #-----Create 'predictions', 'probabilities' lists for descriptor set
    
    # TODO at the moment, some of the columns are factors and need to be converted
    # to numeric
    all.preds <- as.data.frame(apply(all.preds, 2,
                                     function(x) as.numeric(as.character(x))))
    # TODO make the lists data frames if they only have one element
    if (classify) {
      rownames(all.probs) <- IDS
      des.probs.ls <- c(des.probs.ls, list(all.probs))
    }
    
    rownames(all.preds) <- IDS
    des.preds.ls <- c(des.preds.ls, list(all.preds))
    # print(head(all.probs))
    # print(head(all.preds))
  }
  
  names(des.preds.ls) <- cml$des.names
  if (classify)
    names(des.probs.ls) <- cml$des.names
  print(des.preds.ls)
  print(des.probs.ls)
}


# pred.tree <- function( model.data, work.data, classify ){
#   #-----Recursive partitioning using "tree" with splitting criterion deviance and default settings. Specifically:
#   # mincut=5, minimum leaf size
#   # minsize=10, minimum parent size
#   # mindev=.01, within-node deviance must be at least this times that of the root node for node to split
#   cat( "Tree ---------------------------------------", "\n" )
#   n.work<-nrow(work.data)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y") 
#   {
#     work.meth <- tree( as.factor(y)~.,data=model.data, method="class" ) 
#     work.pred <- as.numeric(levels(as.factor(model.data$y)))[ predict( work.meth, work.data, type="class" ) ]
#     work.prob <- predict( work.meth, work.data )[,1]
#   }
#   else 
#   {
#     work.meth <- tree( y~.,data=work.data ) 
#     work.pred <- predict( work.meth, work.data, type="vector" ) 
#   }
#   return( list(pred=work.pred, prob=1-work.prob) )
# }
# 
# #--- RPart has an internal function named "pred.rpart", our function is pred.rpart.chm
# pred.rpart.cml <- function( model.data, work.data, classify ){
#   #-----Recursive partitioning using "rpart" with splitting criterion "information" and
#   # minbucket=5, minimum leaf size
#   # minsplit=10, minimum parent size
#   # maxcompete=0, don't get information on competitive splits
#   # maxsurrogate=0, don't get information on surrogate splits
#   # Possible modifications that have NOT been pursued here:
#   # many ...
#   cat( "RPart --------------------------------------", "\n" )
#   n.work<-nrow(work.data)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y") 
#   {
#     work.meth <- rpart( as.factor(y)~.,data=model.data, method="class", parms=list(split="information"), minsplit=10, minbucket=5, maxcompete=0, maxsurrogate=0 ) 
#     work.pred <- as.numeric(levels(as.factor(model.data$y)))[predict( work.meth, work.data, type="class" ) ]
#     work.prob <- predict( work.meth, work.data )[,1]
#   }
#   else 
#   {
#     work.meth <- rpart( y~.,data=model.data, method="anova", parms=list(split="information"), minsplit=10, minbucket=5, maxcompete=0, maxsurrogate=0 ) 
#     work.pred <- predict( work.meth, work.data, type="vector" ) 
#   }
#   return( list(pred=work.pred, prob=1-work.prob) )
# }
# 
# 
# pred.rf <- function( model.data, work.data, classify){
#   #-----Random Forest using
#   # ntree=100
#   # mtry = sqrt(p) [for classification] and = p/3 [for regression]
#   # nodesize = 5
#   # importance = TRUE (to calculate important descriptors)
#   # Possible modifications that have NOT been pursued here:
#   # many ...
#   cat( "Forest -------------------------------------", "\n" )
#   n.work<-nrow(work.data)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y") 
#   {
#     work.meth <- randomForest( y=as.factor(model.data$y), x=model.data[,-1], ntree=100, nodesize=5 ) 
#     work.pred <- as.numeric(levels(as.factor(model.data$y)))[ predict( work.meth, work.data[,-1] ) ]
#     work.prob <- predict( work.meth, work.data[,-1], type="prob" )[,1]
#   }
#   else 
#   {
#     work.meth <- randomForest( y=model.data$y, x=model.data[,-1], ntree=100, nodesize=5 ) 
#     work.pred <- predict( work.meth, work.data[,-1] ) 
#   }
#   return( list(pred=work.pred, prob=1-work.prob) )
# }
# 
# 
# pred.svm <- function( model.data, work.data, classify ){
#   #-----Support Vector Machine using
#   # kernel = radial.basis
#   # gamma = 1
#   # Possible modifications that have NOT been pursued here:
#   # many ...
#   cat( "SVM ----------------------------------------", "\n" )
#   n.work<-nrow(work.data)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y") 
#   {
#     work.meth <- svm( y=as.factor(model.data$y), x=model.data[,-1], gamma=1, probability=TRUE ) 
#     work.pred <- as.numeric(levels(as.factor(model.data$y)))[ predict( work.meth, work.data[,-1] ) ]
#     temp.prob <- predict( work.meth, work.data[,-1], probability=TRUE )
#     work.prob <- attr(temp.prob, "probabilities")[,1]
#   }
#   else 
#   {
#     work.meth <- svm( y=model.data$y, x=model.data[,-1], gamma=1 ) 
#     work.pred <- predict( work.meth, work.data[,-1] ) 
#   }
#   return( list(pred=work.pred, prob=1-work.prob) )
# }
# 
# 
# pred.nnet <- function( model.data, work.data, classify ){
#   #-----Neural Network using "nnet" with:
#   # size=2
#   # trace=FALSE (don't print out convergence info)
#   # Possible modifications that have NOT been pursued here:
#   # many
#   cat( "NNet ---------------------------------------", "\n" )
#   n.work<-nrow(work.data)
#   work.pred<-rep(NA,n.work)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y") {
#     work.meth <- nnet( as.factor(y)~., data=model.data, size=2, trace=FALSE ) 
#     work.pred <- predict( work.meth, work.data, type="class" )
#     work.prob <- 1 - predict( work.meth, work.data, type="raw" )
#   }
#   return( list(pred=work.pred, prob=1-work.prob) )
# }
# 
# 
# pred.knn <- function( model.data, work.data, classify ){
#   #-----K nearest neighbor using
#   # k = 3
#   # Possible modifications that have NOT been pursued here:
#   # cv option ...
#   cat( "KNN ----------------------------------------", "\n" )
#   n.work<-nrow(work.data)
#   work.pred<-rep(NA,n.work)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y") 
#   {
#     work.meth <- knn( cl=as.factor(model.data$y), train=model.data[,-1], test=work.data[,-1], k=10, prob=TRUE ) 
#     work.pred <- as.numeric(levels(as.factor(model.data$y)))[ work.meth ]
#     work.prob <- abs( work.pred + attr(work.meth, "prob") - 1 )
#   }
#   return( list(pred=work.pred, prob=work.prob) )
# }
# 
# 
# pred.knnflex <- function( model.data, work.data, classify ){
#   #-----K nearest neighbor using
#   # k = 3
#   cat( "KNN ----------------------------------------", "\n" )
#   work.meth<-knn.dist(rbind(model.data[,-1],work.data[,-1]))
#   n.model<-nrow(model.data)
#   n.work<-nrow(work.data)
#   work.prob<-rep(NA,n.work)
#   if (classify=="Y")
#   {
#     work.pred <- as.numeric( knn.predict(train=1:n.model,test=(n.model+1):(n.work+n.model),y=work.data$y,k=10,dist.matrix=work.meth,agg.meth="majority") )
#     work.prob <- knn.predict(train=1:n.model,test=(n.model+1):(n.work+n.model),y=work.data$y,k=10,dist.matrix=work.meth,agg.meth="mean")
#   }
#   else
#   {
#     work.pred <- knn.predict(train=1:n.model,test=(n.model+1):(n.work+n.model),y=work.data$y,k=10,dist.matrix=work.meth,agg.meth="mean")
#   }
#   return( list(pred=work.pred, prob=work.prob) )
# }
# 
# 
# pred.lars <- function( model.data, work.data, classify ){
#   #-----Lars
#   cat( "LARs ---------------------------------------", "\n" )
#   work.meth <- lars( y=model.data[,1],x=as.matrix(model.data[,-1]),type="lar" ) 
#   temp.pred <- predict( work.meth,as.matrix(work.data[,-1]),type="fit" )
#   work.pred <- temp.pred$fit[,dim(temp.pred$fit)[2]]
#   return( list(pred=work.pred) )
# }
# 
# 
# #--------------------------------------------------------------------------------
# pred.ridge <- function( model.data, work.data, classify ){
#   #-----Ridge Regression using
#   # lambda = 0.1
#   cat( "Ridge --------------------------------------", "\n" )
#   work.meth <- lm.ridge( y~.,data=model.data,lambda=0.1 ) 
#   work.pred <- (as.matrix(work.data[,-1])%*%(work.meth$coef/work.meth$scales)) +
#     ((work.meth$ym-(work.meth$xm%*%(work.meth$coef/work.meth$scales)))[1,1])
#   return( list(pred=work.pred) )
# }
# 
# 
# #--------------------------------------------------------------------------------
# pred.enet <- function( model.data, work.data, classify ){
#   #-----enet
#   cat( "ENet ---------------------------------------", "\n" )
#   work.meth <- enet( y=model.data$y, x=as.matrix(model.data[,-1]), lambda=1 ) 
#   temp.pred <- predict( work.meth, as.matrix(work.data[,-1]), type="fit" )
#   work.pred <- temp.pred$fit[,dim(temp.pred$fit)[2]]
#   return( list(pred=work.pred) )
# }
# 
# 
# #--------------------------------------------------------------------------------  
# pred.pcr <- function( model.data, work.data, classify ){
#   #-----Principal components regression using home-grown code 
#   cat( "PCR ----------------------------------------", "\n" )
#   newX <- model.data[,-1]
#   work.meth <- pcr.ZG( X=model.data[,-1], Y=model.data$y, newX=work.data[,-1] )
#   work.pred <- work.meth$Ypred
#   return( list(pred=work.pred) )
# }
# 
# 
# #--------------------------------------------------------------------------------  
# pred.plsR <- function( model.data, work.data, classify ){
#   #-----Partial least squares using home-grown code based on "kernelpls"
#   cat( "PLS ----------------------------------------", "\n" )
#   work.meth <- kernelpls.new( X=model.data[,-1], Y=model.data$y, ncomp=min(nrow(model.data),(ncol(model.data)-1),100), newX=work.data[,-1] )
#   nLV.ZG <- ZhuGhodsi( work.meth$Yvar )
#   work.pred <- work.meth$Ypred[,1,nLV.ZG] 
#   return( list(pred=work.pred) )
# }
# 
# 
# #--------------------------------------------------------------------------------  
# pred.pls.lda.new <- function( model.data, work.data, classify ){
#   #-----Partial least squares using home-grown code based on "kernelpls"
#   cat( "PLSLDA -------------------------------------", "\n" )
#   work.meth <- simpls.lda.new( X=model.data[,-1], Y=model.data$y, ncomp=min(nrow(model.data),(ncol(model.data)-1),100), newX=work.data[,-1] )
#   nLV.ZG <- max( ZhuGhodsi( work.meth$Yvar ), 1 )
#   work.prob <- work.meth$Yprob[,2,nLV.ZG] 
#   work.pred <- work.meth$Ypred[,nLV.ZG] 
#   return( list(pred=work.pred, prob=work.prob) )
# }
# 
