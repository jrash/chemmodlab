#' Fit predictive models to sets of descriptors.
#'
#' \code{ModelTrain} is a generic S3 function that fits a series of 
#' classification or regression
#' models to sets of descriptors and computes cross-validated measures
#' of model performance.
#'
#' @param x a list of numeric descriptor set matrices.  At the moment, only
#' binary and continuous descriptors are supported.  Binary descriptors should
#' be numeric (0 or 1).
#' @param y a numeric vector containing the binary or continuous response.
#' @param d a data frame containing an (optional) ID column,
#' a response column, and descriptor columns.  The columns should be
#' provide in this order.
#' @param ids a logical.  Is an ID column provided?
#' @param xcol.lengths a vector of integers.  It is assumed that the columns
#' in \code{d} are grouped by descriptor set.  The integers specify the
#' number of descriptors in each descriptor set.  They should be ordered as
#' the descriptor sets are ordered in \code{d}.
#' Users can specify multiple descriptor sets. By default there is one
#' descriptor set, namely all columns in \code{d} except the response
#' column and
#' the optional ID column.  Specify \code{xcol.lengths} or \code{xcols},
#' but not both.
#' @param xcols A list of integer vectors.  Each vector contains
#' column indices
#' of \code{data} where a set of descriptor variables is located.
#' Users can specify multiple descriptor sets.  Specify \code{xcol.lengths} or \code{xcols},
#' but not both.
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
#' the parameter values for a model.  The list should have the format of
#' the list constructed by  \code{\link{MakeModelDefaults}}. One can construct
#' a list of parameters using  \code{\link{MakeModelDefaults}} and then
#' modify the parameters.
#' @param verbose verbose mode or not?
#' @param descriptors descriptor sets to compute
#' @param mols molecule file created by rcdk
#' @param ... Additional parameters.
#'
#' @return A list is returned of class \code{\link{chemmodlab}} containing:
#'  \item{all.preds}{a list of lists of data frames.  The elements of the outer
#'   list correspond to each CV split performed by \code{\link{ModelTrain}}. The
#'   elements of the inner list correspond to each descriptor set.  For each
#'   descriptor set and CV split combination, the output is a dataframe
#'   containing all model predictions.  The first column of each data frame
#'   contains the true value of the response.  The remaining columns contain
#'   the predictions for each model.}
#' \item{all.probs}{a list of lists of data frames. Constructed only if there is
#'   a binary response.  The structure is the same as \code{all.preds}, except
#'   that predictions are replaced by "predicted probabilities" (i.e. estimated
#'   probabilities of a response
#'   value of one).  Predicted
#'   probabilities are only reported for classification models.}
#' \item{model.acc}{a list of lists of model accuracy measures.  The elements of
#'   the outer list correspond to each CV split performed by \code{ModelTrain}.
#'   The elements of the inner list correspond to each descriptor set.  For each
#'   descriptor set and CV split combination, a limited collection of 
#'   performance measures are given for each model fit
#'   to the data.  Regression models are assessed with Pearson's \eqn{r} and
#'   \eqn{RMSE}. Classification models are assessed with contingency tables.
#'   For additional model performance measures, see \code{\link{Performance}}}.
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
#' is performed for the specified regression and/or classification
#' models.
#' 
#' Not all modeling strategies will be appropriate for all response
#' types. For example, partial least squares linear discriminant analysis
#' ("PLSLDA")
#' is not directly appropriate for continuous response assays such as
#' percent inhibition, but it can be applied once a threshold value for
#' percent inhibition is used to create a binary (active/inactive) response.
#'
#' See \url{https://jrash.github.io/chemmodlab/} for more 
#' information about the
#' models available (including model default parameters).
#' The default value for argument models includes only some of 
#' the possible values.
#'
#' Sensible default values are selected for each
#' tunable model parameter, however users may set any parameter
#' manually using \code{\link{MakeModelDefaults}} and \code{user.params}.
#'
#' \code{\link{ModelTrain}} predictions are based on k-fold cross-validation,
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
#' @author Jacqueline Hughes-Oliver, Jeremy Ash, Atina Brooks
#' @seealso \code{\link{chemmodlab}}, \code{\link{plot.chemmodlab}},
#'   \code{\link{CombineSplits}},
#'
#' @examples
#' 
#' \dontrun{
#' # A data set with  binary response and multiple descriptor sets
#' data(aid364)
#' 
#' cml <- ModelTrain(aid364, ids = TRUE, xcol.lengths = c(24, 147),
#'                   des.names = c("BurdenNumbers", "Pharmacophores"))
#' cml
#' }
#' 
#' # A continuous response
#' cml <- ModelTrain(USArrests, nsplits = 2, nfolds = 2,
#'                   models = c("KNN", "Lasso", "Tree"))
#' cml
#' 
#' @import methods
#' @import stats
#' @import utils
#' @import rcdk
#' @import fingerprint
#' 
#' @export
ModelTrain <- function(...) UseMethod("ModelTrain")

#' @describeIn ModelTrain Default S3 method
#' @export
ModelTrain.default <- function(x, y, 
                               nfolds = 10,
                               nsplits = 3,
                               seed.in = NA,
                               des.names = NA,
                               models = c("NNet", "PLS", "LAR", "Lasso",
                                          "PLSLDA", "Tree", "SVM", "KNN", "RF"),
                               user.params = NULL, verbose = FALSE,
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
  if (!is(y, "numeric")) stop("'y' must be a numeric vector")

  if (verbose) {
    BackModelTrain(x = x, y = y, xcol.lengths = xcol.lengths,
               nfolds = nfolds, nsplits = nsplits, seed.in = seed.in, 
               des.names = des.names, models = models, user.params = user.params,
               s3method = s3method, verbose = verbose)
  } else {
    suppressWarnings(BackModelTrain(x = x, y = y, xcol.lengths = xcol.lengths,
             nfolds = nfolds, nsplits = nsplits, seed.in = seed.in, 
             des.names = des.names, models = models, user.params = user.params,
             s3method = s3method, verbose = verbose))
  }
}

#' @describeIn ModelTrain S3 method for class 'character'
#' @export
ModelTrain.character <- function(descriptors, y, mols, 
                                 nfolds = 10,
                                 nsplits = 3,
                                 seed.in = NA,
                                 des.names = NA,
                                 models = c("NNet", "PLS", "LAR", "Lasso",
                                            "PLSLDA", "Tree", "SVM", "KNN", "RF"),
                                 user.params = NULL, verbose = FALSE,
                                 ...) {
  
  s3method <- "character"
  if (class(mols[[1]]) != "jobjRef") {
    stop("'mols' should be a molecule list produced by the package rcdk")
  }
  
  if (!all(descriptors %in% c("hybrid", "constitutional", "topological",
                         "electronic", "geometrical", "fp.standard",
                         "fp.extended", "fp.graph", "fp.hybridization",
                         "fp.maccs", "fp.estate", "fp.pubchem", "fp.kr",
                         "fp.shortestpath", "fp.signature", "fp.circular"))) {
    stop("'descriptors' should be a character vector containing descriptors computed by chemmodlab")
  }
  
  desc.ls <- list()
  for (i in 1:length(descriptors)) {
    if(descriptors[i] %in% c("hybrid", "constitutional", "topological",
                           "electronic", "geometrical")) {
      desc.ls[[i]] <- eval.desc(mols, get.desc.names(descriptors[i]))
    } else {
      fp.type <- sub("fp.", "", descriptors[i])
      fing <- lapply(mols, get.fingerprint, type = fp.type)
      desc.ls[[i]] <- fp.to.matrix(fing)
    }
  }
  
  df <- matrix(y, ncol = 1)
  xcol.lengths <- c()
  for (i in 1:length(descriptors)) {
    work.data <- as.matrix(desc.ls[[i]])
    rm.desc <- apply(work.data, 2, function(x) all(is.na(x)) == T)
    if (sum(rm.desc) > 0) {
      if (verbose)
        warning(paste("WARNING.....", sum(rm.desc),
                      "descriptors could not be computed and were omitted"))
      work.data <- subset(work.data, select = !rm.desc)
    }
    rm(rm.desc)
    xcol.lengths <- c(xcol.lengths, ncol(work.data))
    df<- cbind(df, work.data)
  }
  df <- as.data.frame(df)
  
  ModelTrain.data.frame(d = df, ids = F, xcol.lengths = xcol.lengths, verbose = verbose) 
}


#' @describeIn ModelTrain S3 method for class 'data.frame'
#' @export
ModelTrain.data.frame <- function(d,
                            ids = FALSE,
                            xcol.lengths = ifelse(ids,
                                                  length(d) - 2,
                                                  length(d) - 1),
                            xcols = NA,
                            nfolds = 10,
                            nsplits = 3,
                            seed.in = NA,
                            des.names = NA,
                            models = c("NNet", "PLS", "LAR", "Lasso",
                                       "PLSLDA", "Tree", "SVM", "KNN", "RF"),
                            user.params = NULL, verbose = FALSE,
                            ...) {
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
  
  ycol <- ifelse(ids, 2, 1)
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
    # shift up indices by 1 if there is a response column, or 2 if there is a 
    # response and id columns
    for(i in 1:length(xcols)) {
      if(is.na(idcol)) {
        xcols[[i]] <- xcols[[i]] + 1
      } else {
        xcols[[i]] <- xcols[[i]] + 2
      }
    }
  }
  
  # TODO There were still warning messages being produced, so I am just
  # suppressing all warnings when in quiet mode
  if (verbose) {
    BackModelTrain(d = d, ids = ids, xcol.lengths = xcol.lengths, xcols = xcols,
                   nfolds = nfolds, nsplits = nsplits, seed.in = seed.in, 
                   des.names = des.names, models = models, user.params = user.params,
                   s3method = s3method, idcol = idcol, ycol = ycol, verbose = verbose)
  } else {
    suppressWarnings(BackModelTrain(d = d, ids = ids, xcol.lengths = xcol.lengths, xcols = xcols,
                   nfolds = nfolds, nsplits = nsplits, seed.in = seed.in, 
                   des.names = des.names, models = models, user.params = user.params,
                   s3method = s3method, idcol = idcol, ycol = ycol, verbose = verbose))
  }
}

BackModelTrain <- function(x = NA, y = NA, d = NA, ids = NA,
                           xcol.lengths = NA, xcols = NA,
                           idcol = NA, ycol = NA, nfolds, nsplits, seed.in, 
                           des.names, models, user.params, s3method,
                           verbose = FALSE) {
  rn <- rownames(d)
  if (s3method == "data.frame") {
    n.des <- length(xcols)
  } else {
    n.des <- length(x)
  }
  if (!(NA %in% des.names)) {
    if (length(des.names) != n.des)
      stop("'des.names' must be the same length as 'x'")
  } else {
    des.names <- c()
    for (i in 1:n.des) {
      des.names <- c(des.names, paste0("Descriptor Set ", i))
    }
  }
  # checking parameters specified correctly
  if (!(nfolds%%1 == 0) || !is.numeric(nfolds))
    stop("'nfolds' must be an integer")
  # if(length(seed.in) != nsplits)
    # stop("length of 'seed.in' must equal number of splits")
  if (!(nsplits%%1 == 0) || !is.numeric(nsplits))
    stop("'nsplits' must be an integer")
  if (!all(models %in% c("NNet", "PCR", "ENet", "PLS", "Ridge", "LAR", "PLSLDA",
                         "RPart", "Tree", "SVM", "KNN", "RF", "Lasso"))) {
    stop("'models' should be a character vector containing models existing in chemmodlab")
  }

  if (!is.na(seed.in[1])) {
    if (length(seed.in) != nsplits)
      stop("length of 'seed.in' must equal number of splits")
  } else {
    seed.in <- c()
    for (i in 1:nsplits) {
      seed.in <- c(seed.in, as.numeric(paste(rep(i, 5), collapse = "")))
    }
  }

  split.preds.ls <- list()
  split.probs.ls <- list()
  split.model.acc.ls <- list()
  work.data.ls <- list()
  rm.rows <- c()

  for(seed.idx in 1:length(seed.in)){
    
    # all.xcols <- c()
    # for(i in seq_along(xcols)) {
    #   all.xcols <- c(all.xcols, xcols[[i]])
    # }
    # all.xcols <- unique(all.xcols)
    # 
    # all.data <- ReadInData(d, ycol, all.xcols, idcol)[[1]]
    
    des.preds.ls <- list()
    des.probs.ls <- list()
    des.model.acc.ls <- list()
    des.model.obj.ls <- list()

    for (des.idx in 1:n.des) {
      
      current.seed <- seed.in[seed.idx]
      
      # take response and current descriptor set columns
      if (s3method == "data.frame") {
        if (verbose) {
          rid.obj <- ReadInData(d, ycol, xcols[[des.idx]], idcol)
        }
        else {
          suppressWarnings(rid.obj  <- ReadInData(d, ycol, xcols[[des.idx]], idcol))
        }
        work.data <- rid.obj[[1]]
        rm.rows <- c(rm.rows, rid.obj[[3]])
      } else {
        work.data <- data.frame(y = y, x[[des.idx]])
      }
      n.pred <- ncol(work.data) - 1
      n.obs <- nrow(work.data)
    
      #-----Assign folds for nfolds-fold cross-validation
      fold.id <- rep(NA, n.obs)
      num.in.folds <- floor(n.obs/nfolds)
      set.seed(current.seed)
      for (id in 1:(nfolds - 1)) {
        fold.id[sample((1:n.obs)[(is.na(fold.id))], num.in.folds, replace = FALSE)] <- id
      }
      fold.id[(1:n.obs)[(is.na(fold.id))]] <- nfolds
      
      model.acc.ls <- list()

      #-----Determine if we have binary or continuous data
      # if sum == 0 then all of the response is either 0 or 1.  
      if (!exists("classify")) {
        if (sum(!(work.data$y %in% c(1, 0))))
          classify <- F 
        else 
          classify <- T
      }

      #-----Make model parameter list
      if (!is.null(user.params)) {
        params <- MakeModelDefaults(n.obs, n.pred, classify, nfolds)
        params <- SetUserParams(params, user.params)
      } else {
        params <- MakeModelDefaults(n.obs, n.pred, classify, nfolds)
      }

      #-----Start outputing progress to console
      cat("Beginning Analysis for Split: ", seed.idx, "and Descriptor Set: ",
          des.names[des.idx], "\n")
      cat("Number of Descriptors: ", n.pred, "\n")
      cat("Responses: ", n.obs, "\n")
      cat("Starting Seed: ", current.seed, "\n")
      cat("Number of CV folds: ", nfolds, "\n\n")
      all.preds <- data.frame(Observed = work.data$y)
      if (classify)
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
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackTree(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed,  
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...Tree not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackTree(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Tree = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, Tree = work.results$prob)
          model.acc.ls <- c(model.acc.ls, Tree = list(work.results$model.acc))
          model.acc.ls <- c(model.acc.ls, Tree = list(work.results$model.acc))
        }
      }

      #-----'rpart' method
      if (sum(models %in% "RPart") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackRpart(work.data, n.obs,
                                                             n.pred, nfolds, fold.id,
                                                             classify, current.seed,
                                                             params),
                                   error = function(e) {
                                     warning(paste("WARNING...RPart not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackRpart(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, RPart = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, RPart = work.results$prob)
          model.acc.ls <- c(model.acc.ls, RPart = list(work.results$model.acc))
        }
      }

      #-----'randomforest' method
      if (sum(models %in% "RF") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackRf(work.data, n.obs,
                                                          n.pred, nfolds, fold.id,
                                                          classify, current.seed,  
                                                          params),
                                   error = function(e) {
                                     warning(paste("WARNING...RF not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackRf(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, RF = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, RF = work.results$prob)
          model.acc.ls <- c(model.acc.ls, RF = list(work.results$model.acc))
        }
      }

      #-----'svm' method
      if (sum(models %in% "SVM") == 1) {
        work.results <- list()
        pt <- proc.time()
        # TODO you can do this verbose check more simply
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackSvm(work.data, n.obs,
                                                           n.pred, nfolds, fold.id,
                                                           classify, current.seed,  
                                                           params),
                                   error = function(e) {
                                     warning(paste("WARNING...SVM not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackSvm(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, SVM = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, SVM = work.results$prob)
          model.acc.ls <- c(model.acc.ls, SVM = list(work.results$model.acc))
        }
      }

      #-----'nnet' method
      if (sum(models %in% "NNet") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackNnet(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed,  
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...NNet not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackNnet(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, NNet = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, NNet = work.results$prob)
          model.acc.ls <- c(model.acc.ls, NNet = list(work.results$model.acc))
        }
      }

      #-----'knn' method
      if (sum(models %in% "KNN") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackKnn(work.data, n.obs, n.pred,
                                                           nfolds, fold.id, classify,
                                                           current.seed,   
                                                           params),
                                  error = function(e) {
                                    warning(paste("WARNING...KNN not run:", e$message))
                                    work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackKnn(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, KNN = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, KNN = work.results$prob)
          model.acc.ls <- c(model.acc.ls, KNN = list(work.results$model.acc))
        }
      }

      #-----'pls.lda' method
      if ((sum(models %in% "PLSLDA") == 1) && (classify)) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackPlsLdaNew(work.data,
                                                                 n.obs, n.pred, nfolds,
                                                                 fold.id, classify,
                                                                 current.seed),
                                   error = function(e) {
                                     warning(paste("WARNING...PLSLDA not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackPlsLdaNew(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, PLSLDA = work.results$pred)
          if (classify)
            all.probs <- data.frame(all.probs, PLSLDA = work.results$prob)
          model.acc.ls <- c(model.acc.ls, PLSLDA = list(work.results$model.acc))
        }
      }

      #-----'lars' method
      if (sum(models %in% "LAR") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackLars(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify,
                                                            current.seed,   
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...LAR not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackLars(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
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
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackRidge(work.data, n.obs,
                                                             n.pred, nfolds, fold.id,
                                                             classify, current.seed,
                                                             params),
                                   error = function(e) {
                                     warning(paste("WARNING...Ridge not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackRidge(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Ridge = work.results$pred)
          model.acc.ls <- c(model.acc.ls, Ridge = list(work.results$model.acc))
        }
      }
      
      #-----'lm.lasso' method
      if (sum(models %in% "Lasso") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackLasso(work.data, n.obs,
                                                             n.pred, nfolds, fold.id,
                                                             classify, current.seed,
                                                             params),
                                   error = function(e) {
                                     warning(paste("WARNING...Lasso not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackLasso(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Lasso = work.results$pred)
          model.acc.ls <- c(model.acc.ls, Lasso = list(work.results$model.acc))
        }
      }
      
      # TODO suspending lassoGLM until we can get it work right
      # #-----'lm.lassoGLM' method
      # if (sum(models %in% "LassoGLM") == 1) {
      #   work.results <- list()
      #   pt <- proc.time()
      #   st <- system.time(tryCatch(work.results <- BackLassoGLM(work.data, n.obs,
      #                                                        n.pred, nfolds, fold.id,
      #                                                        classify, current.seed,
      #                                                        params),
      #                              error = function(e) {
      #                                warning(paste("WARNING...LassoGLM not run:", e$message))
      #                                work.results <- list()
      #                              }))
      #   PrintTime(pt, st)
      #   if (length(work.results) > 0) {
      #     all.preds <- data.frame(all.preds, LassoGLM = work.results$pred)
      #     model.acc.ls <- c(model.acc.ls, LassoGLM = list(work.results$model.acc))
      #   }
      # }

      #-----'enet' method
      if (sum(models %in% "ENet") == 1) {
        work.results <- list()
        pt <- proc.time()
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackEnet(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed,  
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...ENet not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackEnet(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
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
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackPcr(work.data, n.obs,
                                                           n.pred, nfolds, fold.id,
                                                           classify, current.seed,  
                                                           params),
                                   error = function(e) {
                                     warning(paste("WARNING...PCR not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackPcr(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
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
        if (verbose) {
          st <- system.time(tryCatch(work.results <- BackPlsR(work.data, n.obs,
                                                            n.pred, nfolds, fold.id,
                                                            classify, current.seed,  
                                                            params),
                                   error = function(e) {
                                     warning(paste("WARNING...PLS not run:", e$message))
                                     work.results <- list()
                                       }))
        } else {
          st <- system.time(tryCatch(work.results <- BackPlsR(work.data, n.obs,
                                                    n.pred, nfolds, fold.id,
                                                    classify, current.seed,  
                                                    params), error=function(e){}))
        }
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
      des.model.acc.ls <- c(des.model.acc.ls, list(model.acc.ls))
      if (seed.idx == 1) {
        work.data.ls <- c(work.data.ls,
                          list(work.data[, colnames(work.data) != "y"]))
      }
    }
    names(des.preds.ls) <- des.names
    names(des.model.acc.ls) <- des.names
    if (classify)
      names(des.probs.ls) <- des.names
    split.preds.ls <- c(split.preds.ls, list(des.preds.ls))
    split.model.acc.ls <- c(split.model.acc.ls, list(des.model.acc.ls))
    if (classify)
      split.probs.ls <- c(split.probs.ls, list(des.probs.ls))
  }
  
  # Remove missing observations in all descriptor sets
  # print(rm.rows)
  if(length(rm.rows) > 0) {
    rm.rows <- unique(rm.rows)
    rm.rows.rn <- rn[rm.rows]
    for(seed.idx in 1:length(seed.in)){
      for (des.idx in 1:n.des) {
        work.data.rn <- rownames(split.preds.ls[[seed.idx]][[des.idx]])
        split.preds.ls[[seed.idx]][[des.idx]] <- split.preds.ls[[seed.idx]][[des.idx]][!work.data.rn %in% rm.rows.rn, ]
        if (classify)
          split.prob.ls[[seed.idx]][[des.idx]] <- split.prob.ls[[seed.idx]][[des.idx]][!work.data.rn %in% rm.rows.rn, ]
      }
    }
    y <- work.data$y[!work.data.rn %in% rm.rows.rn]
    } else {
      y <- work.data$y
    }

  cml.result <- chemmodlab(split.preds.ls, split.probs.ls, split.model.acc.ls,
                           classify, y, work.data.ls, params, des.names,
                           models, nsplits)

  return(cml.result)
}



