#' Model parameters for ModelTrain
#'
#' Makes a list containing the default parameters for all
#' models implemented in \code{\link{ModelTrain}}.
#'
#' @param n The number of observations in the data.
#' @param p The number of descriptors in the data.
#' @param classify A logical.  Will classification models be used? (is
#' the response binary?)
#' If false, regression models will be assumed.
#' @param nfolds The number of folds used for k-fold cross validation.
#'
#' @return A list whose elements are dataframes containing the
#' default parameter values for models implemented in
#' \code{\link{ModelTrain}}.
#'
#' @details
#' Sensible default values are selected for each
#' tunable model parameter, however users may set any parameter
#' manually by generating a list with this function and assigning
#' the parameters.
#'
#' See \url{https://jrash.github.io/chemmodlab/} for more information about the
#' models available (including model default parameters).
#'
#' @aliases MakeModelDefaults
#' @author Jeremy Ash
#' @seealso \code{\link{ModelTrain}}, \code{\link{chemmodlab}}
#'
#' @examples
#' params <- MakeModelDefaults(n = nrow(USArrests),
#'  p = ncol(USArrests[, -1]), classify = TRUE, nfolds = 10)
#' params$Forest$mtry <- ncol(USArrests[, -1])-1
#' params
#'
#' cml <- ModelTrain(USArrests, models = "RF", nsplits = 3,
#'  user.params = params)
#'
#'
#' @export
MakeModelDefaults <- function(n, p, classify, nfolds){
  params <- list(NNet = data.frame(size = 2, decay = 0),
                 PCR = NULL,
                 ENet = data.frame(lambda = 1),
                 PLS = data.frame(ncomp = min(floor(n/nfolds), p, 100)),
                 Ridge = data.frame(lambda = .1),
                 LAR = data.frame(lambda = .05),
                 Lasso = data.frame(lambda = .2),
                 # LassoGLM = data.frame(lambda = .1),
                 PLSLDA = data.frame(ncomp = min(floor(n/nfolds), p, 100)),
                 RPart = data.frame(cp = .01),
                 Tree = NULL,
                 SVM = data.frame(gamma = 1, cost = 1, epsilon = .01),
                 KNN = data.frame(k = 10),
                 Forest = data.frame(mtry = if (classify) max(floor(p/3), 1)
                                     else floor(sqrt(p)))
  )
  params
}

# can probably remove the code below

PrintModelDefaults <- function(n, p, classify, nfolds) {
  params <- MakeModelDefaults(n, p, classify, nfolds)
  print(params)
}

SetUserParams <- function(params, user.params){
  for (model in names(user.params)) {
    if (!is.null(user.params[[model]]) && nrow(user.params[[model]]) > 1) {
      stop(cat("Can only set on parameter value for", model))
    }
    grid.set <-
      names(user.params[[model]])[names(user.params[[model]]) %in% names(params[[model]])]
    grid.excess <-
      names(user.params[[model]])[!names(user.params[[model]]) %in% names(params[[model]])]
    for (hpar in grid.set) {
      params[[model]][, hpar] <- user.params[[model]][, hpar]
    }
    if (length(grid.excess) > 0){
      warning(paste0("Unused parameter values for model ", model, ": ", grid.excess))
    }
  }
  params
}
