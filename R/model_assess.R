#' ANOVA and multiple comparisons for chemmodlab objects
#' 
#' \code{CombineSplits} assesses observed performance measures
#' across all splits created by \code{\link{ModelTrain}} and conducts
#' statistical tests to determine the best performing descriptor set and
#' model (D-M) combinations.  
#' 
#' \code{CombineSplits}
#' quantifies how sensitive performance measures are to fold
#' assignments (assignments to training and test sets). 
#' Intuitively, this
#' assesses how much a performance measure may change if a slightly
#' different data set is used.
#' 
#' \code{CombineSplits} is a designed study in that 'experimental' 
#' conditions are defined according to two factors: method (D-M combination) 
#' and split (fold assignment).  As such, an analysis
#' of variance to identify significant differences between 
#' performance measures according to factors and levels is appropriate.
#' 
#' The multiple comparisons plot shows the results for tests for signficance
#' in all pairwise differences of D-M mean performance measures.
#' Because there can be many 
#' estimated mean performance measures for an assay, care must be taken 
#' to adjust for
#' multiple testing, and we do this using the Tukey-Kramer multiple
#' comparison procedure (see Tukey (1953) and Kramer (1956)).
#' If you are having trouble viewing all the components of the plot,
#' make the plotting window larger.
#' 
#' By default, \code{CombineSplits} uses initial enhancement
#' proposed by Kearsley et al. (1996) to assess model performance. 
#' Enhancement at \code{n} tests is the hit
#' rate at \code{n} tests (accumulated actives at \code{n} tests
#' divided by \code{n}, the 
#' number of tests) divided by the proportion of actives in the entire 
#' collection. It is a relative measure of hit rate improvement offered
#' by the new method beyond what can be expected under random selection,
#' and values much larger than one are desired. Initial enhancement is
#' typically taken to be enhancement at \code{n}=300 tests.
#' 
#' Root mean squared error (\code{RMSE}), despite its popularity
#' in statistics, may be  
#' inappropriate for continuous chemical assay responses because
#' it assumes losses 
#' are equal for both under-predicting and over-predicting biological 
#' activity.  A suitable alternative may be initial \code{enhancement}.
#' Other options are the coeffcient of determination (\code{R2})
#' and Spearman's \code{rho}.
#' 
#' For binary chemical assay responses, alternatives to 
#' misclassification rate (\code{error rate}) 
#' (which may be inappropriate because it assigns equal weights to false
#' positives and false negatives) include \code{sensitivity},
#' \code{specificity},
#' area under the receiver operating characteristic curve (\code{ROC}),
#' and initial \code{enhancement}.
#' 
#' @aliases CombineSplits
#' @author Jacqueline Hughes-Oliver, Jeremy Ash, Atina Brooks
#' @seealso \code{\link{chemmodlab}}, \code{\link{ModelTrain}}
#' 
#' @param cml.result an object of class \code{\link{chemmodlab}}.
#' @param at the number of tests to use for initial \code{enhancement}. If 
#' \code{at} is not specified, 
#'  use \code{floor(min(300,n/4))}, where n is the number of compounds.
#' @param metric the model performance measure to use.  This should be
#' one of \code{"error rate"}, \code{"enhancement"}, \code{"R2"},
#' \code{"rho"}, \code{"auc"}, \code{"sensitivity"}, \code{"specificity"}.
#' @param thresh if the predicted probability that a binary response is 
#' 1 is above this threshold, an observation is classified as 1. Used
#' to compute \code{"error rate"}, \code{"sensitivity"}, and 
#' \code{"specificity"}
#' 
#' @references 
#' Kearsley, S.K., Sallamack, S., Fluder, E.M., Andose, J.D., Mosley, R.T.,
#' and Sheridan, R.P. (1996). Chemical similarity using physiochemical
#' property descriptors, J. Chem. Inf. Comput. Sci. 36, 118-127.
#' 
#' Kramer, C. Y. (1956). Extension of multiple range tests to group means
#' with unequal numbers of replications. Biometrics 12, 307-310.
#' 
#' Tukey, J. W. (1953). The problem of multiple comparisons. Unpublished
#' manuscript. In The Collected Works of John W. Tukey VIII. Multiple
#' Comparisons: 1948-1983, Chapman and Hall, New York.
#' 
#' @examples
#' # A data set with  binary response and multiple descriptor sets
#' cml <- ModelTrain(aid364, ids = TRUE, xcol.lengths = c(24, 147), 
#'                   des.names = c("BurdenNumbers", "Pharmacophores"))
#' CombineSplits(cml)
#' 
#' # A continuous response
#' cml <- ModelTrain(USArrests)
#' CombineSplits(cml)
#' 
#' @export

CombineSplits <- function(cml.result, metric = "enhancement",
                          at = NA, thresh = 0.5) {
  
  file <- NA
  y <- cml.result$responses
  # take initial enhancement at 300 tests or if less then 300, the number of y's/4
  if (is.na(at)) {
    at <- min(300, ceiling(length(y)/4))
  } else if (at > length(y)) {
    stop("'at' needs to be smaller than the number of responses")
  }
  
  # makes desciptor set names shorter so that they fit on the MCS plot
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
  if (!cml.result$classify) {
    # then no classification model
    # TODO: Why are they concerned about these objects
    # existing?
    if (exists("sp"))
      rm(sp)
    if (exists("ds"))
      rm(ds)
    if (exists("me"))
      rm(me)
    if (exists("ma"))
      rm(ma)
    
    # TODO Check model accuracy measures
    for (split in 1:length(cml.result$all.preds)) {
      pred <- cml.result$all.preds[[split]]
      for (i in 1:length(pred)) {
        desc <- abbrev.names[i]
        for (j in 2:ncol(pred[[i]])) {
          if (metric == "enhancement") {
            model.acc <- EnhancementCont(pred[[i]][, j], y, at)
          } else if (metric == "R2") {
            model.acc <- (cor(y, pred[[i]][, j]))^2
          } else if (metric == "RMSE") {
            model.acc <- sqrt(sum(y - pred[[i]][, j])^2)
          } else if (metric == "rho") {
            model.acc <- cor(y, pred[[i]][, j], method = "spearman")
          } else {
            stop("y is continuous. 'metric' should be a model accuracy measure
                 implemented for continuous response in ChemModLab")
          }
          colnm <- names(pred[[i]])[j]
          # TO DO do we still need this?
          period <- regexpr(".", colnm, fixed = TRUE)[1]
          if (period > 0)
            meth <- substr(colnm, 1, (period - 1))
          else meth <- colnm
          if (!exists("sp"))
            sp <- split
          else sp <- c(sp, split)
          if (!exists("ds"))
            ds <- desc
          else ds <- c(ds, desc)
          if (!exists("me"))
            me <- meth
          else me <- c(me, meth)
          if (!exists("ma"))
            ma <- model.acc
          else ma <- c(ma, model.acc)
          }
      }
    }
  } else {
    if (exists("sp"))
      rm(sp)
    if (exists("ds"))
      rm(ds)
    if (exists("me"))
      rm(me)
    if (exists("ma"))
      rm(ma)
    for (split in 1:length(cml.result$all.preds)) {
      prob <- cml.result$all.probs[[split]]
      pred <- cml.result$all.preds[[split]]
      for (i in 1:length(prob)) {
        desc <- abbrev.names[i]
        for (j in 2:ncol(prob[[i]])) {
          # TODO: Determine whether to use predicted probabilities or predictions
          # for models that provide both - sometimes they differ.
          yhat <- prob[[i]][, j] > thresh
          if (metric == "enhancement") {
            model.acc <- Enhancement(prob[[i]][, j], y, at)
          } else if (metric == "auc") {
            model.acc <- as.numeric(pROC::auc(y, prob[[i]][, j]))
          } else if (metric == "error rate") {
            model.acc <- mean(y != yhat)
          } else if (metric == "specificity") {
            idx <- y == 0
            model.acc <- mean(y[idx] == yhat[idx])
          } else if (metric == "sensitivity") {
            idx <- y == 1
            model.acc <- mean(y[idx] == yhat[idx])
            # TODO: Remove rho?
            #           } else if (metric == "rho") {
            #             model.acc <- cor(y, prob[[i]][, j], method = "spearman")
          } else {
            stop("y is binary. 'metric' should be a model accuracy measure
                 implemented for binary response in ChemModLab")
          }
          colnm <- names(prob[[i]])[j]
          period <- regexpr(".", colnm, fixed = TRUE)[1]
          if (period > 0)
            meth <- substr(colnm, 1, (period - 1))
          else meth <- colnm
          if (!exists("sp"))
            sp <- split else sp <- c(sp, split)
          if (!exists("ds"))
            ds <- desc else ds <- c(ds, desc)
          if (!exists("me"))
            me <- meth else me <- c(me, meth)
          if (!exists("ma"))
            ma <- model.acc else ma <- c(ma, model.acc)
          }
      }
      for (i in 1:length(pred)) {
        desc <- abbrev.names[i]
        for (j in 2:ncol(pred[[i]])) {
          yhat <- pred[[i]][, j] > thresh
          if (metric == "enhancement") {
            model.acc <- Enhancement(pred[[i]][, j], y, at)
          } else if (metric == "auc") {
            model.acc <- as.numeric(pROC::auc(y, pred[[i]][, j]))
          } else if (metric == "error rate") {
            model.acc <- mean(y != yhat)
          } else if (metric == "specificity") {
            idx <- y == 0
            model.acc <- mean(y[idx] == yhat[idx])
          } else if (metric == "sensitivity") {
            idx <- y == 1
            model.acc <- mean(y[idx] == yhat[idx])
          } else if (metric == "rho") {
            model.acc <- cor(y, pred[[i]][, j], method = "spearman")
          } else {
            stop("y is binary. 'metric' should be a model accuracy measure
                 implemented for binary response in ChemModLab")
          }
          colnm <- names(pred[[i]])[j]
          # TODO: do we still need this?
          period <- regexpr(".", colnm, fixed = TRUE)[1]
          if (period > 0)
            meth <- substr(colnm, 1, (period - 1)) else meth <- colnm
          if (!exists("sp"))
            sp <- split else sp <- c(sp, split)
          if (!exists("ds"))
            ds <- desc else ds <- c(ds, desc)
          if (!exists("me"))
            me <- meth else me <- c(me, meth)
          if (!exists("ma"))
            ma <- model.acc else ma <- c(ma, model.acc)
          }
      }
    }
  }
  
  out <- data.frame(sp, ds, me, ma)
  names(out) <- c("Split", "Descriptor", "Method", "Model.Acc")
  
  # for performance reasons it would be preferable not to calculate the enhancement
  # for the classification methods, but for speed of implementation this is being
  # done here.  Dropping all the model accuracy calculated using predictions from
  # classification methods
  
  out <- out[!duplicated(out[, 1:3]), ]
  out$Split <- factor(out$Split)
  
  methods <- as.character(unique(out$Method))
  num.enhs <- dim(out)[1]
  descriptors <- abbrev.names
  out$Trmt <- vector(mode = "logical", length = num.enhs)
  for (i in 1:num.enhs) {
    out$Trmt[i] <- which(methods == out$Method[i])
  }
  for (i in 1:num.enhs) {
    out$Trmt[i] <- out$Trmt[i] + 100 * (which(descriptors == out$Descriptor[i]))
  }
  
  out$Trmt <- factor(out$Trmt)
  SplitAnova(out, metric, file)
}

SplitAnova <- function(splitdata, metric, file = NA) {
  
  single.desc <- (length(unique(splitdata$Descriptor)) == 1)
  # Calculate anova
  out <- glm(Model.Acc ~ Split + Trmt, data = splitdata)
  # total df - df for residuals
  mod.df <- out$df.null - out$df.residual
  # is this TSS - SSR?
  mod.dev <- out$null.deviance - out$deviance
  mod.ms <- mod.dev/mod.df
  err.df <- out$df.residual
  err.dev <- out$deviance
  err.ms <- err.dev/err.df
  tot.df <- out$df.null
  tot.dev <- out$null.deviance
  f.stat <- mod.ms/err.ms
  p.val <- df(f.stat, mod.df, err.df)
  if (p.val < 1e-04)
    p.val <- "<.0001"
  form.df <- format(c("DF", mod.df, err.df, tot.df), width = 3, justify = "right")
  form.dec.num <- format(c(mod.dev, err.dev, tot.dev, mod.ms, err.ms, f.stat),
                         digits = 4, nsmall = 3)
  form.dec <- format(c("SS", "MS", "F", form.dec.num), digits = 4, nsmall = 3,
                     justify = "right")
  form.p <- format(c("p-value", if (is.numeric(p.val)) round(p.val, digits = 4) else p.val),
                   digits = 1, nsmall = 4, justify = "right")
  cat(paste0("   Analysis of Variance on: '",metric,"'\n"))
  if (single.desc)
    cat(paste("            Using factors: Split and Method\n"))
  else cat(paste(" Using factors: Split and Descriptor/Method combination\n"))
  cat(paste("Source", form.df[1], form.dec[1], form.dec[2], form.dec[3], form.p[1],
            "\n", sep = "   "))
  cat(paste("Model ", form.df[2], form.dec[4], form.dec[7], form.dec[9], form.p[2],
            "\n", sep = "   "))
  cat(paste("Error ", form.df[3], form.dec[5], form.dec[8], "\n", sep = "   "))
  cat(paste("Total ", form.df[4], form.dec[6], "\n", sep = "   "))
  
  # if (!is.na(file))
  #   pdf(file)
  
  # plot(c(-1,1),c(-1,1),xlab='',ylab='',axes=FALSE,type='n')
  # text(-1,.9,str1,adj=c(0,.5),family='mono',cex=.7)
  # text(-1,.8,str2,adj=c(0,.5),family='mono',cex=.7)
  # text(rep(-1,4),seq(.6,.3,by=-.1),c(str3,str4,str5,str6),adj=c(0,0),family='mono',cex=.7)
  root.mse <- sqrt(err.ms)
  all.mean <- mean(splitdata$Model.Acc)
  coef.var <- root.mse/all.mean * 100
  r.sq <- mod.dev/tot.dev
  form.stat.num <- format(c(r.sq, coef.var, root.mse, all.mean), digits = 3, nsmall = 4)
  form.stat <- format(c("R-Square", "Coef Var", "Root MSE", "Mean", form.stat.num),
                      digits = 3, nsmall = 4, justify = "right")
  cat(paste("   ", form.stat[1], form.stat[2], form.stat[3], form.stat[4], "\n",
            sep = "   "))
  cat(paste("   ", form.stat[5], form.stat[6], form.stat[7], form.stat[8], "\n",
            sep = "   "))
  # text(rep(-1,2),seq(.1,0,by=-.1),c(str7,str8),adj=c(0,0),family='mono',cex=.7)
  
  aout <- anova(out)
  f.df <- aout$Df[2]
  f.dev <- aout$Deviance[2]
  f.ms <- f.dev/f.df
  f.f.stat <- f.ms/err.ms
  t.df <- aout$Df[3]
  t.dev <- aout$Deviance[3]
  t.ms <- t.dev/t.df
  t.f.stat <- t.ms/err.ms
  f.p.val <- df(f.f.stat, f.df, t.df)
  if (f.p.val < 1e-04)
    f.p.val <- "<.0001"
  t.p.val <- df(t.f.stat, t.df, err.df)
  if (t.p.val < 1e-04)
    t.p.val <- "<.0001"
  form.df <- format(c("DF", f.df, t.df), width = 3, justify = "right")
  form.dec.num <- format(c(f.dev, f.ms, f.f.stat, t.dev, t.ms, t.f.stat), digits = 3,
                         nsmall = 3)
  form.dec <- format(c("SS", "MS", "F", form.dec.num), digits = 3, nsmall = 3,
                     justify = "right")
  form.p <- format(c("p-value", f.p.val, t.p.val), digits = 1, nsmall = 4, justify = "right")
  form.p <- format(c("p-value", if (is.numeric(f.p.val)) round(f.p.val, digits = 4) else f.p.val,
                     if (is.numeric(t.p.val)) round(t.p.val, digits = 4)
                     else t.p.val), digits = 1, nsmall = 4, justify = "right")
  cat(paste("Source   ", form.df[1], form.dec[1], form.dec[2], form.dec[3], form.p[1],
            "\n", sep = "   "))
  cat(paste("Split    ", form.df[2], form.dec[4], form.dec[5], form.dec[6], form.p[2],
            "\n", sep = "   "))
  if (single.desc)
    cat(paste("Method   ", form.df[3], form.dec[7], form.dec[8], form.dec[9],
              form.p[3], "\n", sep = "   "))
  else cat(paste("Desc/Meth", form.df[3], form.dec[7], form.dec[8], form.dec[9],
                 form.p[3], "\n", sep = "   "))
  # text(rep(-1,3),seq(-.2,-.4,by=-.1),c(str9,str10,str11),adj=c(0,0),family='mono',cex=.7)
  
  # Calculate lsmeans
  lsmeans <- c()
  fit <- fitted.values(out)
  for (i in levels(splitdata$Trmt)) lsmeans <- c(lsmeans, mean(fit[splitdata$Trmt ==
                                                                     i]))
  
  # Calculate Tukey-Kramer adjusted p-values for diff means
  aov.out <- aov(Model.Acc ~ Split + Trmt, data = splitdata)
  tout <- TukeyHSD(aov.out, "Trmt")
  
  # Create matrix of p-values
  n <- length(levels(splitdata$Trmt))
  pval <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      pval[i, j] <- tout$Trmt[n * (i - 1) - i * (i - 1)/2 + j - i, 4]
      pval[j, i] <- pval[i, j]
    }
  }
  
  # Plot
  McsPlot(lsmeans, pval, as.numeric(levels(splitdata$Trmt)),
          as.character(unique(splitdata$Descriptor)),
          as.character(unique(splitdata$Method)), single.desc, metric)
  if (!is.na(file))
    dev.off()
}


