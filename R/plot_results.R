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
#' accumulates \code{y} so that \code{\eqn{\sum_{i=1}^n y_i}{\sum y_i}}
#' is the sum of the
#' \code{y}
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
#' @param cml.result an object of class \code{\link{chemmodlab}}.
#' @param max.select the maximum number of tests to plot for the
#'  accumulation curve. If max.select not specified, use \code{floor(min(300,n/4))},
#'  where n is the number of compounds.
#' @param splits a numeric vector containing the indices of the splits to use to construct
#' accumulation curves.  Default is to use all splits. \code{NA} means the first series
#' of plots are not generate. See \code{Details}.
#' @param models a character vector with models implemented in \code{chemmodlab}.  The 
#' models to use for the second series plots.  See \code{Details}
#' 
#' @author Jacqueline Hughes-Oliver, Jeremy Ash
#' @seealso \code{\link{chemmodlab}}, \code{\link{ModelTrain}}
#' @references Modified from code originally written by
#'   William J. Welch 2001-2002
#'   
#' @examples
#' # A data set with  binary response and multiple descriptor sets
#' 
#' cml <- ModelTrain(aid364, ids = T, xcol.lengths = c(24, 147), 
#'                   des.names = c("BurdenNumbers", "Pharmacophores"))
#' plot(cml)
#' 
#' # A continuous response
#' 
#' cml <- ModelTrain(USArrests)
#' plot(cml)
#' @export
plot.chemmodlab <- function(cml.result, max.select = NA, splits = 1:cml.result$nsplits,
                            models = cml.result$models) {
  
  if (!all(models %in% c("NNet", "PCR", "ENet", "PLS", "Ridge", "LARs", "PLSLDA",
                         "RPart", "Tree", "SVM", "KNN", "Forest", "Forest70", "TreeEns",
                         "RPartEns",
                         "KNNEns"))) {
    stop("'models' should be a character vector containing models existing in chemmodlab")
  }
  
  if (missing(max.select))
    max.select <- min(300,(length(cml.result$responses)/4))

  nsplit <- length(cml.result$all.preds)
  
  # makes desciptor set names shorter so that they fit on the plots
  abbrev.names <- c()
  num.desc <- length(cml.result$des.names)
  des.names <- cml.result$des.names
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
  
  if (!NA %in% splits) {
    for (splidx in splits) {
      preds <- cml.result$all.preds[[splidx]]
      titles <- paste0("Split ", splidx, " : ", gsub("_", " ", names(preds)))
  
      y <- cml.result$responses
      if (sum(!(y %in% c(1, 0))))
        classify <- "N" else classify <- "Y"
      if (classify == "Y")
        num.actives <- sum(y)
      if (classify == "Y")
        probs <- cml.result$all.probs[[splidx]]
  
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
    }
           
    # Plot all descriptors for each method
    if (!NA %in% models) {
      meths <- models
      # for (i in length(titles)) meths <- c(meths, names(preds[[i]])[-1])
      meths <- unique(meths)
      if (classify == "Y") {
        pmeths <- c()
        for (i in length(titles)) pmeths <- c(pmeths, names(probs[[i]])[-1])
        pmeths <- unique(pmeths)
        pmeths <- meths[meths %in% pmeths]
        meths <- meths[!(meths %in% pmeths)]
        for (j in 1:length(pmeths)) {
          p <- data.frame()
          p.labels <- c()
          for (i in 1:length(titles)) {
            if (sum(names(probs[[i]]) %in% pmeths[j]) == 1) {
              if (ncol(p) == 0) {
                p <- data.frame(probs[[i]][, pmeths[j]])
              } else {
                p <- data.frame(p, probs[[i]][, pmeths[j]])
              }
              names(p)[ncol(p)] <- abbrev.names[i]
              p.labels <- c(p.labels, abbrev.names[i])
            }
          }
          if (ncol(p)>0) {
            HitCurve(p, y=y, title=paste("Split", splidx, ":", pmeths[j]),
                     phat.labels=p.labels)
          }
        }
      }
      for (j in 1:length(meths)) {
        p <- data.frame()
        p.labels <- c()
        for (i in 1:length(titles)) {
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

#' @export
plotSingleRun <- function(cml.result, splitidx, desidx) {
  y <- cml.result$response
  desc <- cml.result$data[[desidx]]

  # is this necessary, wasnt data cleaning already done?
  y <- y[apply(desc, 1, function(x) sum(is.na(x)) == 0)]
  desc <- subset(desc, apply(desc, 1, function(x) sum(is.na(x)) == 0))
  desc <- subset(desc, !is.na(y))
  y <- y[!is.na(y)]
  if (sum(!(y %in% c(1, 0))))
    classify <- "N" else classify <- "Y"
  preds <- cml.result$all.preds[[splitidx]][[desidx]]
  if (classify == "Y")
    probs <- cml.result$all.probs[[splitidx]][[desidx]]
  if ((classify == "Y") && (ncol(probs) > 1)) {
    HitCurve(probs[, -1], y = probs[, 1], phat.labels = names(probs)[-1])
    ContCurve(preds[, !(names(preds) %in% names(probs))],
              yhat.labels = names(preds)[!(names(preds) %in%
                                             names(probs))], y = preds[, 1],
              curves.only = TRUE, start.col = (ncol(probs) - 1))
  } else {
    ContCurve(preds[, -1], y = preds[, 1], yhat.labels = names(preds)[-1])
  }
  if ((classify == "Y") && (ncol(probs) > 1)) {
    HitCurve(probs[, -1], y = probs[, 1], phat.labels = names(probs)[-1])
    ContCurve(preds[, !(names(preds) %in% names(probs))],
              yhat.labels = names(preds)[!(names(preds) %in%
                                             names(probs))], y = preds[, 1],
              curves.only = TRUE, start.col = (ncol(probs) - 1))
  } else {
    ContCurve(preds[, -1], y = preds[, 1], yhat.labels = names(preds)[-1])
  }
}
