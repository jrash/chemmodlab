#' Plot method for the chemmodlab class.
#'
#' \code{plot.chemmodlab} takes a \code{\link{chemmodlab}} object output by the
#' \code{\link{ModelTrain}} function and creates a series of accumulation curve
#' plots for assesing model and descriptor set performance.
#'
#' @details
#' For a binary response, the accumulation curve plots the number of assay hits
#' identified as a function of the number of tests conducted, where testing
#' order is determined by the predicted probability of a response being positive
#' obtained from k-fold cross
#' validation. Given a particular compound collection, larger accumulations are
#' preferable.
#'
#' The accumulation curve has also been extended to continuous responses.
#' Assuming large positive values of a continuous response y are preferable,
#' ChemModLab
#' accumulates \code{y} so that \eqn{\sum y_i} is the sum of the \code{y}
#' over the first \code{n} tests. This extension includes the binary-response
#' accumulation curve as a special case.
#'
#' By default, we display accumulation curves up to 300 tests, not for the
#' entire collection, to focus on the goal of finding actives as early as
#' possible.
#'
#' There are two main series of plots generated:
#'
#' @section First plot series:
#'  There is one plot per CV split and descriptor set
#'  combination. The accumulation curve for each model is compared.
#'
#' @section Second plot series:
#'  There is one plot per CV split and model fit. The
#'  accumulation curve for each descriptor set is
#'  compared.
#'
#' @aliases plot.chemmodlab
#' @param x an object of class \code{\link{chemmodlab}}.
#' @param max.select the maximum number of tests to plot for the
#'  accumulation curve. If \code{max.select} is not specified, 
#'  use \code{floor(min(300,n/4))},
#'  where n is the number of compounds.
#' @param splits a numeric vector containing the indices of the splits to use to construct
#' accumulation curves.  Default is to use all splits. \code{NA} means the first series
#' of plots are not generated. See \code{Details}.
#' @param meths a character vector with statistical methods implemented in
#' \code{chemmodlab}.  The 
#' statistical methods to use for the second series of plots.  See \code{Details}.
#' @param ... other parameters to be passed through to plotting functions.
#' 
#' @author Jacqueline Hughes-Oliver, Jeremy Ash, Atina Brooks
#' @seealso \code{\link{chemmodlab}}, \code{\link{ModelTrain}}
#' @references Modified from code originally written by
#'   William J. Welch 2001-2002
#'   
#' @examples
#' # A data set with  binary response and multiple descriptor sets
#' cml <- ModelTrain(aid364, ids = TRUE, xcol.lengths = c(24, 147), 
#'                   des.names = c("BurdenNumbers", "Pharmacophores"))
#' plot(cml)
#' 
#' # A continuous response
#' cml <- ModelTrain(USArrests)
#' plot(cml)
#' @export
plot.chemmodlab <- function(x, max.select = NA, splits = 1:x$nsplits,
                            meths = x$models, ...) {
  
  # This function will modify graphical parameters
  # Reset old parameters upon exiting
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  if (!all(meths %in% c("NNet", "PCR", "ENet", "PLS", "Ridge", "LAR", "PLSLDA",
                         "RPart", "Tree", "SVM", "KNN", "RF"))) {
    stop("'meths' should be a character vector containing methods existing in chemmodlab")
  }
  
  if (missing(max.select))
    max.select <- min(300,(length(x$responses)/4))

  nsplit <- length(x$all.preds)
  
  # Makes desciptor set names shorter so that they fit on the plots
  abbrev.names <- c()
  num.desc <- length(x$des.names)
  des.names <- x$des.names
  if (grepl("Descriptor Set", des.names[1])) {
    for (i in 1:num.desc) {
      # TO DO add option to specify abbreviated names?
      abbrev.names <- c(abbrev.names, paste0("Des", i))
    }
  } else {
    for (i in 1:num.desc) {
      abbrev.names <- c(abbrev.names, substr(des.names[i], 1, 4))
    }
  }
  
  for (splidx in splits) {
    preds <- x$all.preds[[splidx]]
    titles <- paste0("Split ", splidx, " : ", gsub("_", " ", names(preds)))

    y <- x$responses
    if (sum(!(y %in% c(1, 0)))) classify <- "N" else classify <- "Y"
    if (classify == "Y")
      num.actives <- sum(y)
    if (classify == "Y")
      probs <- x$all.probs[[splidx]]

    # TO DO why removing KNN, etc predictions?
    if (classify == "N") {
      for (i in length(titles))
        preds[[i]] <- preds[[i]][, !(names(preds[[i]]) %in% c("KNN", "NNet", "PLSLDA"))]
    }
    
    #Plot all methods for each descriptor set
    for (i in 1:num.desc) {
      if ((classify == "Y") && (ncol(probs[[i]]) > 1)) {
        HitCurve(probs[[i]][, -1], y = y, title = titles[i],
                 phat.labels = names(probs[[i]])[-1])
        ContCurve(preds[[i]][, !(names(preds[[i]]) %in% names(probs[[i]]))],
                  y = y, curves.only = TRUE, start.col = (ncol(probs[[i]]) - 1),
                  title = titles[i],
                  yhat.labels =
                    names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))])
      } else {
        ContCurve(preds[[i]][, -1], y = y, title = titles[i],
                  yhat.labels = names(preds[[i]])[-1])
      }
    }
    
    # Plot all descriptor accumulation curves for each method
    # DONT make these plots if there is only one descriptor set
    if (num.desc > 1){
      # Only use the models that were succesfully fit to the data
      # these will be the models with columns in the all.preds dataframes
      pred.meths <- c()
      for (i in length(titles)) pred.meths <- c(pred.meths, names(preds[[i]])[-1])
      meths <- meths[(meths %in% pred.meths)]
      
      # Make sure method names are unique
      meths <- unique(meths)
      
      if (classify == "Y") {
        prob.meths <- c()
        for (i in length(titles)) prob.meths <- c(prob.meths, names(probs[[i]])[-1])
        prob.meths <- unique(prob.meths)
        prob.meths <- meths[meths %in% prob.meths]
        meths <- meths[!(meths %in% prob.meths)]
        for (j in 1:length(prob.meths)) {
          p <- data.frame()
          p.labels <- c()
          for (i in 1:length(titles)) {
            if (sum(names(probs[[i]]) %in% prob.meths[j]) == 1) {
              if (ncol(p) == 0) {
                p <- data.frame(probs[[i]][, prob.meths[j]])
              } else {
                p <- data.frame(p, probs[[i]][, prob.meths[j]])
              }
              names(p)[ncol(p)] <- abbrev.names[i]
              p.labels <- c(p.labels, abbrev.names[i])
            }
          }
          if (ncol(p)>0) {
            HitCurve(p, y=y, title=paste("Split", splidx, ":", prob.meths[j]),
                     phat.labels=p.labels)
          }
        }
      }
      
      for (j in seq_along(meths)) {
        p <- data.frame()
        p.labels <- c()
        for (i in seq_along(titles)) {
          if (sum(names(preds[[i]]) %in% meths[j]) == 1) {
            if (ncol(p) == 0) {
              p <- data.frame(preds[[i]][, meths[j]])
            } else {
              p <- data.frame(p, preds[[i]][, meths[j]])
            }
            names(p)[ncol(p)] <- abbrev.names[i]
            p.labels <- c(p.labels, abbrev.names[i])
          }
        }
        if (ncol(p) > 0) {
          ContCurve(p, y = y, title = paste("Split", splidx, ":", meths[j]),
                    yhat.labels = p.labels)
        }
      }
    }
  }
}
