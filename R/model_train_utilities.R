#' Constructor for the chemmodlab object
#'
#' Constructor for the chemmodlab object
#' @param all.preds a list of lists of dataframes.  The elements of the outer
#'   list correspond to each split performed by \code{\link{ModelTrain}}. The
#'   elements of the inner list correspond to each descriptor set.  For each
#'   descriptor set and CV split combination, the output is a dataframe
#'   containing all model predictions.  The first column of each data frame
#'   contains the true value of the response.  The remaining columns correspond
#'   to the models fit to the data.
#' @param all.probs a list of lists of dataframes. Constructed only if there is
#'   a binary response.  The structure is the same as \code{all.preds}, except
#'   that predictions are replaced by predicted probabilities of a response
#'   value of one.  Predicted
#'   probabilities are only reported for classification models (see
#'   \code{\link{ModelTrain}})
#' @param model.acc a list of lists of model accuracy measures.  The elements of
#'   the outer list correspond each split performed by \code{\link{ModelTrain}}.
#'   The elements of the inner list correspond to each descriptor set.  For each
#'   descriptor set and CV split combination model accuracy measures for each model fit
#'   to the data.  Regression models are assessed with Pearson's \eqn{r} and
#'   \eqn{RMSE} Classification models are assessed with contingency tables.
#' @param classify a logical.  Was classification models used for binary
#'   response?
#' @param responses a numeric vector.  The true value of the response.
#' @param data a list of numeric matrices.  Each matrix is a descriptor set used
#'   as model input. The first column of each matrix is the response vector, and
#'   the remaining columns are descriptors.
#' @param params a list of dataframes as made by
#'   \code{\link{MakeModelDefaults}}.  Each dataframe contains the parameters to
#'   be set for a particular model.
#'
#' @export

chemmodlab <- function(all.preds, all.probs, model.acc, classify, responses,
                       data, params) {
  cml <- structure(list(all.preds = all.preds, all.probs = all.probs,
                        model.acc = model.acc,
                        classify = classify, responses = responses, data = data,
                        params = params), class = "chemmodlab")
}

# Prints the prediction accuracy on held out data: confusion matrix for
# classification and pearson correlation and Root MSE for regression
BackAssess <- function(obs.resp, pred.resp, imp.desc, classify) {
  if (classify == "Y") {
    return(table(obs.resp, pred.resp))
  } else {
    acc_vec <- c(Pearson_Correlation = cor(obs.resp, pred.resp),
                 Root_MSE = sqrt(sum((obs.resp - pred.resp)^2)))
    return(acc_vec)
  }
}
