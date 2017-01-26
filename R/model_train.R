#' Fit predictive models to sets of descriptors.
#'
#' \code{ModelTrain} fits a series of classification or regression
#' models to sets of descriptors and computes cross-validated measures
#' of model performance.
#'
#' @param data a data frame containing an (optional) ID column,
#' a response column, and descriptor columns.  The columns should be
#' provide in this order.
#' The response variable should either be binary
#' (represented as a numeric vector with 0 or 1 values) or
#' continuous.  At the moment, only numeric descriptors are supported.
#' @param ids a logical.  Is an ID column provided?
#' @param xcol.lengths a vector of integers.  It is assumed that the columns
#'  in \code{data} are grouped by descriptor set.  The integers specify the
#'  number of descriptors in each descriptor set.  They should be ordered as
#'  the descriptor sets are ordered in \code{data}.
#'  Users can specify multiple descriptor sets. By default there is one
#'  descriptor set, all columns in \code{data} except \code{ycol} and
#'  \code{idcol}.
#' @param nfolds the number of folds to use for each cross
#' validation split.
#' @param nsplits the number of splits to use for repeated
#' cross validation.
#' @param seed.in a numeric vector with length equal to \code{nsplits}.
#' The seeds are used to randomly assign folds to observations for each
#' repeated cross-validation split. If \code{NA}, the first seed will be 
#' 11111, the second will be 22222, and so on. 
#' @param models a character vector specifying the regression or
#' classification models to use.  The strings must match models
#' implemented in `chemmodlab` (see Details).
#' @param des.names a character vector specifying the names for each
#' descriptor
#' set.  The length of the vector must match the number of descriptor sets.
#' If \code{NA}, each descriptor set will be named "Descriptor Set i", where
#' i is the number of the descriptor set.
#' @param user.params a list of data frames where each data frame contains
#' the parameters values for a model.  The list should have the format of
#' the list constructed by \code{\link{MakeModelDefaults}}. One can construct
#' a list of parameters using \code{\link{MakeModelDefaults}} and then
#' modify the parameters.
#'
#' @return A list is returned of class \code{\link{chemmodlab}} containing:
#'  \item{all.preds}{a list of lists of dataframes.  The elements of the outer
#'   list correspond to each CV split performed by \code{\link{ModelTrain}}. The
#'   elements of the inner list correspond to each descriptor set.  For each
#'   descriptor set and CV split combination, the output is a dataframe
#'   containing all model predictions.  The first column of each data frame
#'   contains the true value of the response.  The remaining columns contain
#'   the predictions for each model.}
#' \item{all.probs}{a list of lists of dataframes. Constructed only if there is
#'   a binary response.  The structure is the same as \code{all.preds}, except
#'   that predictions are replaced by "predicted probabilities" (ie. estimated
#'   probabilities of a response
#'   value of one).  Predicted
#'   probabilities are only reported for classification models.}
#' \item{model.acc}{a list of lists of model accuracy measures.  The elements of
#'   the outer list correspond each CV split performed by \code{ModelTrain}.
#'   The elements of the inner list correspond to each descriptor set.  For each
#'   descriptor set and CV split combination model accuracy measures for each model fit
#'   to the data.  Regression models are assessed with Pearson's \eqn{r} and
#'   \eqn{RMSE}. Classification models are assessed with contingency tables.}
#' \item{classify}{a logical.  Were classification models used for binary
#'   response?}
#' \item{responses}{a numeric vector.  The observed value of the response.}
#' \item{data}{a list of numeric matrices.  Each matrix is a descriptor set used
#'   as model input.}
#' \item{params}{a list of data frames as made by
#'   \code{\link{MakeModelDefaults}}.  Each data frame contains the parameters to
#'   be set for a particular model.}
#' \item{des.names}{a character vector specifying the descriptor set names.  NA if 
#'   unspecified.}
#' \item{models}{a character vector specifying the models fit to the data.}
#' \item{nsplits}{number of CV splits performed.}
#'
#' @details
#' Multiple descriptor sets can be specified
#' by the user. For each descriptor set, repeated k-fold cross validation
#' is performed for the spcified number of regression and/or classification
#' models.
#' 
#' Not all modeling strategies will be appropriate for all response
#' types. For example, partial least squares linear discriminant analysis
#' ("PLSLDA")
#' is not directly appropriate for continuous response assays such as
#' percent inhibition, but it can be applied once a threshold value for
#' percent inhibition is used to create a binary (active/inactive) response.
#'
#' See \url{https://pages.github.ncsu.edu/jrash/chemmodlab/} for more 
#' information about the
#' models available (including model default parameters).
#'
#' Sensible default values are selected for each
#' tunable model parameter, however users may set any parameter
#' manually using \code{\link{MakeModelDefaults}} and \code{user.params}.
#'
#' \code{ModelTrain} predictions are based on k-fold cross-validation,
#' where the dataset is randomly divided into k parts, each containing
#' approximately equal numbers of compounds. Treating one of these parts
#' as a "test set" the remaining
#' k-1 parts are combined together as a "training set"
#' and used to build a model from the desired modeling technique and
#' descriptor set. This model is then applied to the "test set" to obtain
#' predictions. The process is repeated, holding out each of the k parts
#' in turn. One advantage of k-fold cross-validation is reduction in bias
#' from using the same data to both build and assess a model. Another
#' advantage is the increased precision of error estimation offered by
#' k-fold cross validation over a one-time split.
#'
#' Recognizing that the definition of folds in k-fold cross validation
#' may have an impact on the observed performance measures, all models
#' are built using the same definition of folds. This process is repeated
#' to obtain multiple separate k-fold cross validation runs resulting in
#' multiple separate definitions of folds.  The number of these "splits"
#' is specified by \code{nsplits}.
#'
#' Observed performance measures are
#' assessed across all splits using \code{\link{CombineSplits}}.  This
#' function assesses how sensitive performance measures are to fold
#' assignments, or changes to the training and test sets. 
#' Statistical tests are used to determine the best performing model and
#' descriptor set combination.
#'
#' @aliases ModelTrain
#' @author Jacqueline Hughes-Oliver, Jeremy Ash
#' @seealso \code{\link{chemmodlab}}, \code{\link{plot.chemmodlab}},
#'   \code{\link{CombineSplits}},
#' @references ?
#'
#' @examples
#' 
#' # A data set with  binary response and multiple descriptor sets
#' 
#' cml <- ModelTrain(aid364, ids = T, xcol.lengths = c(24, 147), 
#'                   des.names = c("BurdenNumbers", "Pharmacophores"))
#' cml
#' 
#' # A continuous response
#' 
#' cml <- ModelTrain(USArrests)
#' cml
#'
#' @export
ModelTrain <- function(data,
                       ids = F,
                       xcol.lengths = ifelse(ids,
                                      length(data) - 2,
                                      length(data) - 1),
                       nfolds = 10,
                       nsplits = 3,
                       seed.in = NA,
                       des.names = NA,
                       models = c("NNet", "PLS", "LARs",
                                  "PLSLDA", "Tree", "SVM", "KNN", "Forest"),
                       user.params = NULL) {
  # TO DO: probably can get rid of the idcol parameter

  # checking parameters specified correctly
  if (!is(data, "data.frame"))
    stop("'data' must be a data frame")
  if (!(nfolds%%1 == 0) || !is.numeric(nfolds))
    stop("'nfolds' must be an integer")
  # if(length(seed.in) != nsplits)
    # stop("length of 'seed.in' must equal number of splits")
  if (!(nsplits%%1 == 0) || !is.numeric(nsplits))
    stop("'nsplits' must be an integer")
  if (!is(xcol.lengths, "vector"))
    stop("'xcol.lengths' must be a vector of integers") else {
      for (i in 1:length(xcol.lengths)) {
        if (!(xcol.lengths[[i]]%%1 == 0) || !is.numeric(xcol.lengths[[i]]))
          stop("'xcol.lengths' must be a list of integers")
      }
    }
  if (!all(models %in% c("NNet", "PCR", "ENet", "PLS", "Ridge", "LARs", "PLSLDA",
                         "RPart", "Tree", "SVM", "KNN", "Forest", "Forest70", "TreeEns",
                         "RPartEns",
                         "KNNEns"))) {
    stop("'models' should be a character vector containing models existing in chemmodlab")
  }
  if (!is(ids, "logical")) {
    stop("'ids' should be a logical")
  }
  
  ycol <- ifelse(ids, 2, 1)
  idcol <- ifelse(ids, 1, NA)

  if (!(NA %in% des.names)) {
    if (length(des.names) != length(xcol.lengths))
      stop("'des.names' must be the same length as 'xcol.lengths'")
  } else {
    des.names <- c()
    for (i in 1:length(xcol.lengths)) {
      des.names <- c(des.names, paste0("Descriptor Set ", i))
    }
  }

  if (!is.na(seed.in[1])) {
    if (length(seed.in) != nsplits)
      stop("length of 'seed.in' must equal number of splits")
  } else {
    seed.in <- seq(11111, as.numeric(paste(rep(nsplits, 5), collapse = "")),
                   by = 11111)
  }

  split.preds.ls <- list()
  split.probs.ls <- list()
  split.model.acc.ls <- list()
  work.data.ls <- list()

  n.obs <- nrow(data)

  Funcs <- lsf.str()
  for(seed.idx in 1:length(seed.in)){

    current.seed <- seed.in[seed.idx]

    #-----Assign folds for nfolds-fold cross-validation
    fold.id <- rep(NA, n.obs)
    num.in.folds <- floor(n.obs/nfolds)
    set.seed(current.seed)
    for (id in 1:(nfolds - 1)) {
      fold.id[sample((1:n.obs)[(is.na(fold.id))], num.in.folds, replace = FALSE)] <- id
    }
    fold.id[(1:n.obs)[(is.na(fold.id))]] <- nfolds

    # use the descriptor column numbers to find the corresponding columns in the 
    # data frame
    xcols <- list()
    xcols[[1]] <- 1:xcol.lengths[1]
    for(i in 2:length(xcol.lengths)){
      l1 <- xcols[[i-1]][xcol.lengths[i-1]]
      l2 <- xcol.lengths[i]
      xcols[[i]] <- (l1+1):(l1+l2)
    }
    # shift up indices by 1 if there is a response column, or 2 if there is a 
    # response and id columns
    for(i in 1:length(xcols)){
      if(is.na(idcol)){
        xcols[[i]] <- xcols[[i]] + 1
      } else {
        xcols[[i]] <- xcols[[i]] + 2
      }
    }
    
    des.preds.ls <- list()
    des.probs.ls <- list()
    des.model.acc.ls <- list()

    for (des.idx in 1:length(xcol.lengths)) {
      model.acc.ls <- list()
      work.data <- ReadInData(data, ycol, xcols[[des.idx]], idcol)[[1]]
      n.pred <- ncol(work.data) - 1

      #-----Determine if we have binary or continuous data
      # if sum == 0 then all of the response is either 0 or 1.  This is interpreted as
      # F go to else
      if (!exists("classify")) {
        if (sum(!(work.data$y %in% c(1, 0))))
          classify <- "N" else classify <- "Y"
      }

      #-----Make model parameter list

      if (!is.null(user.params)) {
        params <- MakeModelDefaults(n.obs, n.pred, ifelse(classify == "Y", T, F), nfolds)
        params <- SetUserParams(params, user.params)
      } else {
        params <- MakeModelDefaults(n.obs, n.pred, ifelse(classify == "Y", T, F), nfolds)
      }

      #-----Start outputing progress to console
      cat("Begining Analysis for Split: ", seed.idx, "and Descriptor Set: ",
          des.names[des.idx], "\n")
      cat("Number of Descriptors: ", n.pred, "\n")
      cat("Responses: ", n.obs, "\n")
      cat("Starting Seed: ", current.seed, "\n")
      cat("Number of CV folds: ", nfolds, "\n\n")
      all.preds <- data.frame(Observed = work.data$y)
      if (classify == "Y")
        all.probs <- data.frame(Observed = work.data$y)
      else all.probs <- NA
      all.imp.descs <- data.frame(Descriptors = names(work.data))

      IDS <- rownames(work.data)

      #### This section is calling the machine learning functions and saving the results.
      #### If results are succesfully returned, the predictions are appended to the
      #### prediction data frame

      #-----'tree' method
      # dont use method if requested more than once?
      if (sum(models %in% "Tree") == 1) {
        work.results <- list()
        # get how much time the current process has used
        pt <- proc.time()
        # get the time that expression used
        st <- system.time(tryCatch(work.results <- BackTree(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed, nperm,
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...Tree not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Tree = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, Tree = work.results$prob)
          model.acc.ls <- c(model.acc.ls, Tree = list(work.results$model.acc))
        }
      }

      #-----'rpart' method
      if (sum(models %in% "RPart") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackRpart(work.data, n.obs,
                                                             n.pred, nfolds, fold.id,
                                                             classify, current.seed,
                                                             nperm, params),
                                   error = function(e) {
                                     warning(paste("WARNING...RPart not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, RPart = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, RPart = work.results$prob)
          model.acc.ls <- c(model.acc.ls, RPart = list(work.results$model.acc))
        }
      }

      #-----'randomforest' method
      if (sum(models %in% "Forest") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackRf(work.data, n.obs,
                                                          n.pred, nfolds, fold.id,
                                                          classify, current.seed, nperm,
                                                          params),
                                   error = function(e) {
                                     warning(paste("WARNING...RF not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, RF = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, RF = work.results$prob)
          model.acc.ls <- c(model.acc.ls, RF = list(work.results$model.acc))
        }
      }

      #-----'svm' method
      if (sum(models %in% "SVM") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackSvm(work.data, n.obs,
                                                           n.pred, nfolds, fold.id,
                                                           classify, current.seed, nperm,
                                                           params),
                                   error = function(e) {
                                     warning(paste("WARNING...SVM not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, SVM = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, SVM = work.results$prob)
          model.acc.ls <- c(model.acc.ls, SVM = list(work.results$model.acc))
        }
      }

      #-----'nnet' method
      if ((sum(models %in% "NNet") == 1) && (classify == "Y")) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackNnet(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed, nperm,
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...NNet not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, NNet = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, NNet = work.results$prob)
          model.acc.ls <- c(model.acc.ls, NNet = list(work.results$model.acc))
        }
      }

      #-----'knn' method
      if ((sum(models %in% "KNN") == 1) && (classify == "Y")) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackKnn(work.data, n.obs, n.pred,
                                                           nfolds, fold.id, classify,
                                                           current.seed, nperm, params),
                                  error = function(e) {
                                    warning(paste("WARNING...KNN not run:", e$message))
                                    work.results <- list()
                                  }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, KNN = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, KNN = work.results$prob)
          model.acc.ls <- c(model.acc.ls, KNN = list(work.results$model.acc))
        }
      }

      #-----'pls.lda' method
      if ((sum(models %in% "PLSLDA") == 1) && (classify == "Y")) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackPlsLdaNew(work.data,
                                                                 n.obs, n.pred, nfolds,
                                                                 fold.id, classify,
                                                                 current.seed, nperm),
                                   error = function(e) {
                                     warning(paste("WARNING...PLSLDA not run:", e$message))
                                     work.results <- list()
                                   }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, PLSLDA = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, PLSLDA = work.results$prob)
          model.acc.ls <- c(model.acc.ls, PLSLDA = list(work.results$model.acc))
        }
      }

      #-----'lars' method
      if (sum(models %in% "LARs") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackLars(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify,
                                                            current.seed, nperm, params),
                                   error = function(e) {
                                     warning(paste("WARNING...LAR not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, LAR = work.results$pred)
          model.acc.ls <- c(model.acc.ls, LAR = list(work.results$model.acc))
        }
      }

      #-----'lm.ridge' method
      if (sum(models %in% "Ridge") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackRidge(work.data, n.obs,
                                                             n.pred, nfolds, fold.id,
                                                             classify, current.seed,
                                                             nperm, params),
                                   error = function(e) {
                                     warning(paste("WARNING...Ridge not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Ridge = work.results$pred)
          model.acc.ls <- c(model.acc.ls, Ridge = list(work.results$model.acc))
        }
      }

      #-----'enet' method
      if (sum(models %in% "ENet") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackEnet(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed, nperm,
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...ENet not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, ENet = work.results$pred)
          model.acc.ls <- c(model.acc.ls, ENet = list(work.results$model.acc))
        }
      }

      #-----'PcrZG' method
      if (sum(models %in% "PCR") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackPcr(work.data, n.obs,
                                                           n.pred, nfolds, fold.id,
                                                           classify, current.seed, nperm,
                                                           params),
                                   error = function(e) {
                                     warning(paste("WARNING...PCR not run:", e$message))
                                     work.results <- list()
                                     }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, PCR = work.results$pred)
          model.acc.ls <- c(model.acc.ls, PCR = list(work.results$model.acc))
        }
      }

      #-----'pls.R' method
      if (sum(models %in% "PLS") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackPlsR(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed, nperm,
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...PLS not run:", e$message))
                                     work.results <- list()
                                     }))
        # why saving empty list to work.results replace output.to.log with warning()
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, PLS = work.results$pred)
          model.acc.ls <- c(model.acc.ls, PLS = list(work.results$model.acc))
        }
      }

      cat("Ending Analysis for Split: ", seed.idx, "and Descriptor Set: ",
          des.names[des.idx], "\n\n\n")

      #-----Create 'predictions', 'probabilities' and 'model accuracy' lists 
      # for descriptor set

      # TO DO at the moment, some of the columns are factors and need to be converted
      # to numeric
      all.preds <- as.data.frame(apply(all.preds, 2,
                                       function(x) as.numeric(as.character(x))))
      # TO DO make the lists data frames if they only have one element
      if (classify == "Y") {
        rownames(all.probs) <- IDS
        des.probs.ls <- c(des.probs.ls, list(all.probs))
      }

      rownames(all.preds) <- IDS
      des.preds.ls <- c(des.preds.ls, list(all.preds))
      des.model.acc.ls <- c(des.model.acc.ls, list(model.acc.ls))
      if (seed.idx == 1)
        work.data.ls <- c(work.data.ls,
                          list(work.data[, colnames(work.data) != "y"]))
    }
    names(des.preds.ls) <- des.names
    names(des.model.acc.ls) <- des.names
    if (classify == "Y")
      names(des.probs.ls) <- des.names
    split.preds.ls <- c(split.preds.ls, list(des.preds.ls))
    split.model.acc.ls <- c(split.model.acc.ls, list(des.model.acc.ls))
    if (classify == "Y")
      split.probs.ls <- c(split.probs.ls, list(des.probs.ls))
  }

  cml.result <- chemmodlab(split.preds.ls, split.probs.ls, split.model.acc.ls,
                           classify, work.data$y, work.data.ls, params, des.names,
                           models, nsplits)

  return(cml.result)
}
