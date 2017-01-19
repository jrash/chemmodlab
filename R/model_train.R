#' Fit predictive models to sets of descriptors.
#'
#' \code{ModelTrain} fits a series of classification or regression
#' models to sets of descriptors and computes cross-validated measures
#' of model performance.
#'
#' @param data A data frame containing a response column, descriptor columns
#' and an (optional) id column.
#' @param ycol The index of the response column.
#' @param xcols A list of integer vectors.  Each vector contains
#'  column indices
#'  of \code{data} where a set of descriptor variables is located.
#'  Users can specify multiple descriptor sets. By default there is one
#'  descriptor set, all columns in \code{data} except \code{ycol} and
#'  \code{idcol}.
#' @param idcol The index of the column of \code{data} containing
#' the compound labels. By default this is \code{NA}, and the row labels
#' of \code{data} will be used.
#' @param nfolds The number of folds to use for each cross
#' validation split.
#' @param nsplits The number of splits to use for repeated
#' cross validation.
#' @param seed.in A numeric vector with length equal to \code{nsplits}.
#' The seeds are used to randomly assign folds to observations for each
#' repeated cross-validation split.
#' @param models A character vector specifying the regression or
#' classification models to use.  The strings must match models
#' implemented in `chemmodlab` (see Details).
#' @param des.names A character vector specifying the names for each
#' descriptor
#' set.  The length of the vector must match the number of descriptor sets.
#' If \code{NA}, descriptor set i will be named "Descriptor Set i".
#' @param user.params A list of data frames where each data frame contains
#' the parameters values for a model.  The list should have the format of
#' the list constructed by \code{\link{MakeModelDefaults}}. One can construct
#' a list of parameters using \code{\link{MakeModelDefaults}} and then
#' modify the parameters.
#'
#' @return A list is returned of class \code{\link{chemmodlab}} containing:
#'  \item{all.preds}{a list of lists of data frames.  The elements of the outer
#'   list correspond to each split performed by ModelTrain. The
#'   elements of the inner list correspond to each descriptor set.  The elements
#'   of the inner list are data frames containing all model predictions for a
#'   split and descriptor set combination.  The first column of each data frame
#'   contains the true value of the response.  The remaining columns contain
#'   the predictions for each model.}
#' \item{all.probs}{a list of lists of data frames. Constructed only if there is
#'   a binary response.  The structure is the same as \code{all.preds}, except
#'   that predictions are replaced by predicted probabilities.  Predicted
#'   probabilities are only reported for classification models.}
#' \item{model.acc}{a list of lists of model accuracy measures.  The elements of
#'   the outer list correspond to each split performed by ModelTrain.
#'   The elements of the inner list correspond to each descriptor set.  The
#'   inner list contains model accuracy measures for each model fit
#'   to the data.  Regression models are assessed with Pearson's \eqn{R^2} and
#'   \eqn{RMSE}. Classification models are assessed with contingency tables.}
#' \item{classify}{a logical.  Were classification models used for binary
#'   response?}
#' \item{responses}{a numeric vector.  The true value of the response.}
#' \item{data}{a list of numeric matrices.  Each matrix is a descriptor set used
#'   as model input. The first column of each matrix is the response vector, and
#'   the remaining columns are descriptors.}
#' \item{params}{a list of data frames as made by
#'   \code{\link{MakeModelDefaults}}.  Each data frame contains the parameters to
#'   be set for a particular model.}
#'
#'
#' @details
#' \code{ModelTrain} fits a number of model fitting operations for each
#' data set provided. The response variable should either be binary or
#' continuous.  Descriptors can be any object as long as the models used
#' can take the object as input.
#'
#' Not all modeling strategies will be appropriate for all response
#' types. For example, partial least squares linear discriminant analysis
#' ("PLSLDA")
#' is not directly appropriate for continuous response assays such as
#' percent inhibition, but it can be applied once a threshold value for
#' percent inhibition is used to create a binary (active/inactive) response.
#'
#' See \url{https://jrash.github.io/ChemModLab/} for more information about the
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
#' Multiple descriptor sets can be specified
#' by the user. For each descriptor set, repeated k-fold cross validation
#' is performed for the spcified number of regression and/or classification
#' models.
#'
#' @aliases ModelTrain
#' @author Jacqueline Hughes-Oliver, Jeremy Ash
#' @seealso \code{\link{chemmodlab}}, \code{\link{plot.chemmodlab}},
#'   \code{\link{CombineSplits}},
#' @references ?
#' @examples
#' cml <- ModelTrain(USArrests, nsplits = 3)
#' cml
#'
#' bin <- rbinom(50, 1, .1)
#' cml <- ModelTrain(cbind(bin, USArrests), nsplits = 3)
#' cml
#' @export
ModelTrain <- function(data,
                       ycol = ifelse(is.na(idcol), 1, 2),
                       xcols = ifelse(is.na(idcol),
                                      list(seq(2, ncol(data))),
                                      list(seq(3, ncol(data)))),
                       idcol = NA,
                       nfolds = 10,
                       nsplits = 1,
                       seed.in = NA,
                       des.names = NA,
                       models = c("NNet", "PLS", "LARs",
                                  "PLSLDA", "Tree", "SVM", "KNN", "Forest"),
                       user.params = NULL) {
  # TO DO: probably can get rid of the idcol parameter

  #-----Background benchmarking program
  # yfilein=input response data: response in column ycol xfilein=input descriptors
  # data: descriptors in columns xcols

  # checking parameters specified correctly

  if (!is(data, "data.frame"))
    stop("'data' must be a list of data frames or matrices")
  if (!(ycol%%1 == 0))
    stop("'ycol' must be an integer")
  if (!(nfolds%%1 == 0))
    stop("'nfolds' must be an integer")
  # if(length(seed.in) != nsplits)
    # stop("length of 'seed.in' must equal number of splits")
  if (!(nsplits%%1 == 0))
    stop("'nsplits' must be an integer")
  if (!is(xcols, "list"))
    stop("'xcols' must be a list of integers") else {
      for (i in 1:length(xcols)) {
        if (!(all.equal(xcols[[i]], as.integer(xcols[[i]]))))
          stop("'xcols' must be a list of integers")
      }
    }
  if (!all(models %in% c("NNet", "PCR", "ENet", "PLS", "Ridge", "LARs", "PLSLDA",
                         "RPart", "Tree", "SVM", "KNN", "Forest", "Forest70", "TreeEns", "RPartEns",
                         "KNNEns"))) {
    stop("'models' should be a character vector containing models existing in chemmodlab")
  }
  if (!is.na(idcol) && !(idcol%%1 == 0)) {
    stop("'idcol' should be an integer or NA")
  }
  meths <- models

  if (!is.na(des.names)) {
    if (length(des.names) != length(xcols))
      stop("'des.names' must be the same length as 'xcols'")
  } else {
    des.names <- c()
    for (i in 1:length(xcols)) {
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

    des.preds.ls <- list()
    des.probs.ls <- list()
    des.model.acc.ls <- list()

    for (des.idx in 1:length(xcols)) {
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
      if (sum(meths %in% "Tree") == 1) {
        work.results <- list()
        # get how much time the current process has used
        pt <- proc.time()
        # get the time that expression used
        st <- system.time(tryCatch(work.results <- BackTree(work.data, n.obs,
                                                            n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                              warning(paste("WARNING...Tree not run:", e$message))
                                                              work.results <- list()
                                                            }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Tree = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, Tree = work.results$prob)
          model.acc.ls <- c(model.acc.ls, Tree = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, Tree=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'rpart' method
      if (sum(meths %in% "RPart") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackRpart(work.data, n.obs,
                                                             n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                               warning(paste("WARNING...RPart not run:", e$message))
                                                               work.results <- list()
                                                             }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, RPart = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, RPart = work.results$prob)
          model.acc.ls <- c(model.acc.ls, RPart = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, RPart=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'randomforest' method
      if (sum(meths %in% "Forest") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackRf(work.data, n.obs,
                                                          n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                            warning(paste("WARNING...RF not run:", e$message))
                                                            work.results <- list()
                                                          }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, RF = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, RF = work.results$prob)
          model.acc.ls <- c(model.acc.ls, RF = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, RF=as.numeric(names(work.data) %in%
          # work.results$impdesc) )
        }
      }

      #-----'svm' method
      if (sum(meths %in% "SVM") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackSvm(work.data, n.obs,
                                                           n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                             warning(paste("WARNING...SVM not run:", e$message))
                                                             work.results <- list()
                                                           }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, SVM = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, SVM = work.results$prob)
          model.acc.ls <- c(model.acc.ls, SVM = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, SVM=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'nnet' method
      if ((sum(meths %in% "NNet") == 1) && (classify == "Y")) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackNnet(work.data, n.obs,
                                                            n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                              warning(paste("WARNING...NNet not run:", e$message))
                                                              work.results <- list()
                                                            }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, NNet = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, NNet = work.results$prob)
          model.acc.ls <- c(model.acc.ls, NNet = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, NNet=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'knn' method
      if ((sum(meths %in% "KNN") == 1) && (classify == "Y")) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackKnn(work.data, n.obs,
                                                           n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                             warning(paste("WARNING...KNN not run:", e$message))
                                                             work.results <- list()
                                                           }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, KNN = work.results$pred)
          if (classify == "Y")
            all.probs <- data.frame(all.probs, KNN = work.results$prob)
          model.acc.ls <- c(model.acc.ls, KNN = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, KNN=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'pls.lda' method
      if ((sum(meths %in% "PLSLDA") == 1) && (classify == "Y")) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackPlsLdaNew(work.data,
                                                                 n.obs, n.pred, nfolds, fold.id, classify, current.seed, nperm),
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
          # all.imp.descs <- data.frame( all.imp.descs, PLSLDA=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'lars' method
      if (sum(meths %in% "LARs") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackLars(work.data, n.obs,
                                                            n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                              warning(paste("WARNING...LAR not run:", e$message))
                                                              work.results <- list()
                                                            }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, LAR = work.results$pred)
          model.acc.ls <- c(model.acc.ls, LAR = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, LAR=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'lm.ridge' method
      if (sum(meths %in% "Ridge") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackRidge(work.data, n.obs,
                                                             n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                               warning(paste("WARNING...Ridge not run:", e$message))
                                                               work.results <- list()
                                                             }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, Ridge = work.results$pred)
          model.acc.ls <- c(model.acc.ls, Ridge = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, Ridge=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'enet' method
      if (sum(meths %in% "ENet") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackEnet(work.data, n.obs,
                                                            n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                              warning(paste("WARNING...ENet not run:", e$message))
                                                              work.results <- list()
                                                            }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, ENet = work.results$pred)
          model.acc.ls <- c(model.acc.ls, ENet = list(work.results$model.acc))
          # all.imp.descs <- data.frame( all.imp.descs, ENet=as.numeric(names(work.data) %in% work.results$impdesc) )
        }
      }

      #-----'PcrZG' method
      if (sum(meths %in% "PCR") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackPcr(work.data, n.obs,
                                                           n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                             warning(paste("WARNING...PCR not run:", e$message))
                                                             work.results <- list()
                                                           }))
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, PCR = work.results$pred)
          model.acc.ls <- c(model.acc.ls, PCR = list(work.results$model.acc))
          #  all.imp.descs <- data.frame( all.imp.descs, PCR=as.numeric(names(work.data          #  %in% work.results$impdesc) )
        }
      }

      #-----'pls.R' method
      if (sum(meths %in% "PLS") == 1) {
        work.results <- list()
        pt <- proc.time()
        st <- system.time(tryCatch(work.results <- BackPlsR(work.data, n.obs,
                                                            n.pred, nfolds, fold.id, classify, current.seed, nperm, params), error = function(e) {
                                                              warning(paste("WARNING...PLS not run:", e$message))
                                                              work.results <- list()
                                                            }))
        # why saving empty list to work.results replace output.to.log with warning()
        PrintTime(pt, st)
        if (length(work.results) > 0) {
          all.preds <- data.frame(all.preds, PLS = work.results$pred)
          model.acc.ls <- c(model.acc.ls, PLS = list(work.results$model.acc))
          #  all.imp.descs <- data.frame( all.imp.descs, PLS=as.numeric(names(work.data          #  %in% work.results$impdesc) )
        }
      }

      cat("Ending Analysis for Split: ", seed.idx, "and Descriptor Set: ",
          des.names[des.idx], "\n\n\n")

      #-----Create 'predictions', 'probabilities' and 'model accuracy' lists for descriptor set

      # TO DO at the moment, some of the columns are factors and need to be converted
      # to numeric
      all.preds <- as.data.frame(apply(all.preds, 2, function(x) as.numeric(as.character(x))))
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
                           classify, work.data$y, work.data.ls, params)

  return(cml.result)
}
