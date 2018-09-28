# Need to go back and check the probabilities are reported correctly
PredTree <- function( y, model.data, pred.data, classify, params ){
  #-----Recursive partitioning using "tree" with splitting criterion deviance and default settings. Specifically:
  # mincut=5, minimum leaf size
  # minsize=10, minimum parent size
  # mindev=.01, within-node deviance must be at least this times that of the root node for node to split
  n.work<-nrow(pred.data)
  work.prob<-rep(NA,n.work)
  if (classify) {
    work.meth <- tree::tree( as.factor(y)~.,data=model.data, method="class" ) 
    work.pred <- as.numeric(levels(as.factor(y)))[predict(work.meth, pred.data, type="class")]
    work.prob <- predict(work.meth, pred.data)[, 1]
  }
  else {
    work.meth <- tree::tree( y~.,data=model.data ) 
    work.pred <- predict(work.meth, pred.data, type="vector") 
  }
  return(list(pred=work.pred, prob=1-work.prob))
}

#--- RPart has an internal function named "pred.rpart", our function is pred.rpart.chm
PredRpart <- function( y, model.data, pred.data, classify, params ){
  #-----Recursive partitioning using "rpart" with splitting criterion "information" and
  # minbucket=5, minimum leaf size
  # minsplit=10, minimum parent size
  # maxcompete=0, don't get information on competitive splits
  # maxsurrogate=0, don't get information on surrogate splits
  # Possible modifications that have NOT been pursued here:
  # many ...
  n.work<-nrow(pred.data)
  work.prob<-rep(NA,n.work)
  if (classify) {
    work.meth <- rpart::rpart(as.factor(y) ~ ., data = model.data,
                              method = "class", parms = list(split = "information"),
                              control = rpart::rpart.control(minsplit = 10,
                                                             minbucket = 5, maxcompete = 0,
                                                             maxsurrogate = 0,
                                                             cp = params$RPart$cp))
    work.pred <- as.numeric(levels(as.factor(y)))[predict(work.meth,
                                                          pred.data,
                                                          type = "class")]
    work.prob <- predict(work.meth, pred.data)[, 1]
  }
  else {
    work.meth <- rpart::rpart(y ~ ., data = model.data,
                              method = "anova",
                              control = rpart::rpart.control(minsplit = 10,
                                                             minbucket = 5,
                                                             maxcompete = 0,
                                                             maxsurrogate = 0)) 
    work.pred <- predict(work.meth, pred.data, type="vector") 
  }
  return(list(pred=work.pred, prob=1-work.prob))
}


PredRf <- function( y, model.data, pred.data, classify){
  #-----Random Forest using
  # ntree=100
  # mtry = sqrt(p) [for classification] and = p/3 [for regression]
  # nodesize = 5
  # importance = TRUE (to calculate important descriptors)
  # Possible modifications that have NOT been pursued here:
  # many ...
  cat( "Forest -------------------------------------", "\n" )
  n.work<-nrow(pred.data)
  work.prob<-rep(NA,n.work)
  if (classify) 
  {
    work.meth <- randomForest( y=as.factor(model.data$y), x=model.data[,-1], ntree=100, nodesize=5 ) 
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[ predict( work.meth, pred.data[,-1] ) ]
    work.prob <- predict( work.meth, pred.data[,-1], type="prob" )[,1]
  }
  else 
  {
    work.meth <- randomForest( y=model.data$y, x=model.data[,-1], ntree=100, nodesize=5 ) 
    work.pred <- predict( work.meth, pred.data[,-1] ) 
  }
  return( list(pred=work.pred, prob=1-work.prob) )
}


PredSvm <- function( y, model.data, pred.data, classify, params ){
  #-----Support Vector Machine using
  # kernel = radial.basis
  # gamma = 1
  # Possible modifications that have NOT been pursued here:
  # many ...
  cat( "SVM ----------------------------------------", "\n" )
  n.work<-nrow(pred.data)
  work.prob<-rep(NA,n.work)
  if (classify) 
  {
    work.meth <- svm( y=as.factor(model.data$y), x=model.data[,-1], gamma=1, probability=TRUE ) 
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[ predict( work.meth, pred.data[,-1] ) ]
    temp.prob <- predict( work.meth, pred.data[,-1], probability=TRUE )
    work.prob <- attr(temp.prob, "probabilities")[,1]
  }
  else 
  {
    work.meth <- svm( y=model.data$y, x=model.data[,-1], gamma=1 ) 
    work.pred <- predict( work.meth, pred.data[,-1] ) 
  }
  return( list(pred=work.pred, prob=1-work.prob) )
}


PredNnet <- function( y, model.data, pred.data, classify, params ){
  #-----Neural Network using "nnet" with:
  # size=2
  # trace=FALSE (don't print out convergence info)
  # Possible modifications that have NOT been pursued here:
  # many
  cat( "NNet ---------------------------------------", "\n" )
  n.work<-nrow(pred.data)
  work.pred<-rep(NA,n.work)
  work.prob<-rep(NA,n.work)
  if (classify) {
    work.meth <- nnet( as.factor(y)~., data=model.data, size=2, trace=FALSE ) 
    work.pred <- predict( work.meth, pred.data, type="class" )
    work.prob <- 1 - predict( work.meth, pred.data, type="raw" )
  }
  return( list(pred=work.pred, prob=1-work.prob) )
}


PredKnn <- function( y, model.data, pred.data, classify, params ){
  #-----K nearest neighbor using
  # k = 3
  # Possible modifications that have NOT been pursued here:
  # cv option ...
  cat( "KNN ----------------------------------------", "\n" )
  n.work<-nrow(pred.data)
  work.pred<-rep(NA,n.work)
  work.prob<-rep(NA,n.work)
  if (classify) 
  {
    work.meth <- knn( cl=as.factor(model.data$y), train=model.data[,-1], test=pred.data[,-1], k=10, prob=TRUE ) 
    work.pred <- as.numeric(levels(as.factor(model.data$y)))[ work.meth ]
    work.prob <- abs( work.pred + attr(work.meth, "prob") - 1 )
  }
  return( list(pred=work.pred, prob=work.prob) )
}


PredLars <- function( y, model.data, pred.data, classify, params ){
  #-----Lars
  cat( "LARs ---------------------------------------", "\n" )
  work.meth <- lars( y=model.data[,1],x=as.matrix(model.data[,-1]),type="lar" ) 
  temp.pred <- predict( work.meth,as.matrix(pred.data[,-1]),type="fit" )
  work.pred <- temp.pred$fit[,dim(temp.pred$fit)[2]]
  return( list(pred=work.pred) )
}


#--------------------------------------------------------------------------------
PredRidge <- function( y, model.data, pred.data, classify, params ){
  #-----Ridge Regression using
  # lambda = 0.1
  cat( "Ridge --------------------------------------", "\n" )
  work.meth <- lm.ridge( y~.,data=model.data,lambda=0.1 ) 
  work.pred <- (as.matrix(pred.data[,-1])%*%(work.meth$coef/work.meth$scales)) +
    ((work.meth$ym-(work.meth$xm%*%(work.meth$coef/work.meth$scales)))[1,1])
  return( list(pred=work.pred) )
}


#--------------------------------------------------------------------------------
PredEnet <- function( y, model.data, pred.data, classify, params ){
  #-----enet
  cat( "ENet ---------------------------------------", "\n" )
  work.meth <- enet( y=model.data$y, x=as.matrix(model.data[,-1]), lambda=1 ) 
  temp.pred <- predict( work.meth, as.matrix(pred.data[,-1]), type="fit" )
  work.pred <- temp.pred$fit[,dim(temp.pred$fit)[2]]
  return( list(pred=work.pred) )
}


#--------------------------------------------------------------------------------  
PredPcr <- function( y, model.data, pred.data, classify, params ){
  #-----Principal components regression using home-grown code 
  cat( "PCR ----------------------------------------", "\n" )
  newX <- model.data[,-1]
  work.meth <- pcr.ZG( X=model.data[,-1], Y=model.data$y, newX=pred.data[,-1] )
  work.pred <- work.meth$Ypred
  return( list(pred=work.pred) )
}


#--------------------------------------------------------------------------------  
PredPlsR <- function( y, model.data, pred.data, classify, params ){
  #-----Partial least squares using home-grown code based on "kernelpls"
  cat( "PLS ----------------------------------------", "\n" )
  work.meth <- kernelpls.new( X=model.data[,-1], Y=model.data$y, ncomp=min(nrow(model.data),(ncol(model.data)-1),100), newX=pred.data[,-1] )
  nLV.ZG <- ZhuGhodsi( work.meth$Yvar )
  work.pred <- work.meth$Ypred[,1,nLV.ZG] 
  return( list(pred=work.pred) )
}


#--------------------------------------------------------------------------------  
PredPlsLdaNew <- function( y, model.data, pred.data, classify, params ){
  #-----Partial least squares using home-grown code based on "kernelpls"
  cat( "PLSLDA -------------------------------------", "\n" )
  work.meth <- simpls.lda.new( X=model.data[,-1], Y=model.data$y, ncomp=min(nrow(model.data),(ncol(model.data)-1),100), newX=pred.data[,-1] )
  nLV.ZG <- max( ZhuGhodsi( work.meth$Yvar ), 1 )
  work.prob <- work.meth$Yprob[,2,nLV.ZG] 
  work.pred <- work.meth$Ypred[,nLV.ZG] 
  return( list(pred=work.pred, prob=work.prob) )
}

