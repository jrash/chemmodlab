
#' @export
MakeModelDefaults <- function(n, p, classify, nfolds){
  params <- list(NNet = data.frame(size = 2, decay = 0),
                 PCR = NULL,
                 ENet = data.frame(lambda = 1),
                 PLS = data.frame(ncomp = min(floor(n/nfolds), p, 100)),
                 Ridge = data.frame(lambda = .1),
                 LARs = NULL,
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

#' @export
PrintModelDefaults <- function(n, p, classify, nfolds) {
  params <- MakeModelDefaults(n, p, classify, nfolds)
  print(params)
}

#' @export
SetUserParams <- function(params, user.params){
  for (model in names(user.params)) {
    if (!is.null(user.params[[model]]) && nrow(user.params[[model]]) > 1) {
      stop(cat("Can only set on parameter value for", model))
    }
    grid.set <-
      names(user.params[[model]])[names(user.params[[model]]) %in% names(params[[model]])]
    grid.excess <-
      names(user.params[[model]])[! names(user.params[[model]]) %in% names(params[[model]])]
    for (hpar in grid.set) {
      params[[model]][, hpar] <- user.params[[model]][, hpar]
    }
    if (length(grid.excess) > 0){
      warning(paste0("Unused parameter values for model ", model, ": ", grid.excess))
    }
  }
  params
}
