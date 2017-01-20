#' Plot method for the chemmodlab class.
#'
#' \code{plot.chemmodlab} takes a \code{\link{chemmodlab}} object output by the
#' \code{\link{ModelTrain}} function and creates a series of accumulation curve
#' plots for assesing model and descriptor set performance.
#'
#' For a binary response, the accumulation curve plots the number of assay hits
#' identified as a function of the number of tests conducted, where testing
#' order is determined by the predicted probabilities obtained from k-fold cross
#' validation. Given a particular compound collection, larger accumulations are
#' preferable.
#'
#' The accumulation curve has also been extended to continuous responses.
#' Assuming large values of a continuous response y are preferable, ChemModLab
#' accumulates y so that \eqn{\sum_{i=1}^n y_i}{\sum y_i} is the sum of the y
#' over the first n tests. This extension includes the binary-response
#' accumulation curve as a special case.
#'
#' By default, we display accumulation curves up to 300 tests, not for the
#' entire collection, to focus on the goal of finding actives as early as
#' possible.
#'
#' There are two main series of plots generated:
#'
#' @section First plot series:
#'  There is one plot per split and descriptor set
#'  combination. The accumulation curve for each model is compared.
#'
#' @section Second plot series:
#'  There is one plot per split and model fit. The
#'  accumulation curve for each descriptor set combination is
#'  compared.
#'
#' @aliases plot.chemmodlab
#' @param cml.result an object of class \code{\link{chemmodlab}}.
#' @param max.select the maximum number of tests to plot for the
#'  accumulation curve. If max.select not specified, use \code{floor(min(300,n/4))},
#'  where n is the number of compounds.
#' @author Jacqueline Hughes-Oliver, Jeremy Ash
#' @seealso \code{\link{chemmodlab}}, \code{\link{modelTrain}}, \code{\link{modelTrain}}
#' @references Modified from code originally written by
#'   William J. Welch 2001-2002
#' @examples
#' cml <- modelTrain(USArrests, nsplits = 3)
#' plot(cml)
#'
#' bin <- rbinom(50, 1, .1)
#' cml <- ModelTrain(cbind(bin, USArrests), nsplits = 3)
#' plot(cml)
#'
#' @export
plot.chemmodlab <- function(cml.result, max.select = NA) {
  # bndir,phdir,apdir,fpdir,chdir: contain the dir name of the various descriptor
  #                                set runs, e.g. Split1/BurdenNumbers/2007-02-26-16-18-07,
  #                                note no begining or ending slashes, any of these
  #                                may be missing, at least one required
  # response: the name of the response file, e.g Outcome.csv
  # dest: the dest directory, where the combined files will reside, no
  #       begining or ending slashes, e.g Split1
  # root: the root dir, ends in a slash

  # TO DO: Are we still using this parameter? - only used for hit.maps
  # num.actives: used for continuous responses where the number of actives
  #              can't be easily determined

  #   if (missing("bndir")) bndir=""
  #   if (missing("phdir")) phdir=""
  #   if (missing("apdir")) apdir=""
  #   if (missing("fpdir")) fpdir=""
  #   if (missing("chdir")) chdir=""
  #
  #   i<-1
  #   preds<-list()
  #   if (bndir!="") { preds[[i]]<-read.csv(paste(root,bndir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
  #   i<-i+1 }
  #   if (phdir!="") { preds[[i]]<-read.csv(paste(root,phdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
  #   i<-i+1 }
  #   if (apdir!="") { preds[[i]]<-read.csv(paste(root,apdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
  #   i<-i+1 }
  #   if (fpdir!="") { preds[[i]]<-read.csv(paste(root,fpdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
  #   i<-i+1 }
  #   if (chdir!="") { preds[[i]]<-read.csv(paste(root,chdir,"/pred.csv",sep=""),header=TRUE,skip=1,row.names=1)
  #   i<-i+1 }
  if (missing(max.select))
    max.select <- min(300,(length(cml.result$responses)/4))

  nsplit <- length(cml.result$all.preds)
  for (splidx in 1:nsplit) {
    preds <- cml.result$all.preds[[splidx]]
    titles <- paste0("Split ", splidx, " : ", gsub("_", " ", names(preds)))
    abrev.titles <- c()
    for (j in 1:length(preds)) {
      abrev.titles <- c(abrev.titles, paste0("Des", j))
    }
    #   if (bndir!="") { titles<-c(titles,"Burden Numbers")
    #   abrev.titles<-c(abrev.titles,"B#") }
    #   if (phdir!="") { titles<-c(titles,"Pharmacophores")
    #   abrev.titles<-c(abrev.titles,"Ph") }
    #   if (apdir!="") { titles<-c(titles,"Atom Pairs")
    #   abrev.titles<-c(abrev.titles,"AP") }
    #   if (fpdir!="") { titles<-c(titles,"Fragment Pairs")
    #   abrev.titles<-c(abrev.titles,"FP") }
    #   if (chdir!="") { titles<-c(titles,"Carhart Atom Pairs")
    #   abrev.titles<-c(abrev.titles,"CAP") }

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

    # desc<-data.frame()
    #     if (chdir!="")
    #       desc<-read.csv(paste(root,"Carharts.csv",sep=""),header=TRUE, row.names=1)
    #
    # pdf(paste(root,dest,"/plots.pdf",sep=""))
    #Plot all methods for each descriptor set
    for (i in 1:length(titles)) {
      if ((classify == "Y") && (ncol(probs[[i]]) > 1)) {
        HitCurve(probs[[i]][, -1], y = y, title = titles[i],
                 phat.labels = names(probs[[i]])[-1])
        ContCurve(preds[[i]][, !(names(preds[[i]]) %in% names(probs[[i]]))],
                  y = y, curves.only = TRUE, start.col = (ncol(probs[[i]]) - 1),
                  title = titles[i],
                  yhat.labels =
                    names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))])
        #         if (ncol(desc)==0)
        #           hit.map(data.frame(probs[[i]][,-1],preds[[i]][,!(names(preds[[i]]) %in% names(probs[[i]]))]),
        #                   y=y, title=titles[i], labels=c(names(probs[[i]])[-1],names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))]))
        #         else
        #           hit.map(data.frame(probs[[i]][,-1],preds[[i]][,!(names(preds[[i]]) %in% names(probs[[i]]))]),
        #                   y=y, y.labels=row.names(desc), descriptors=desc, title=titles[i],
        #                   labels=c(names(probs[[i]])[-1],names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))]))

      } else {
        ContCurve(preds[[i]][, -1], y = y, title = titles[i], yhat.labels = names(preds[[i]])[-1])
        #         if (ncol(desc)==0)
        #           hit.map(preds[[i]][,-1], y=y, labels=names(preds[[i]])[-1],
        #                   title=titles[i], max.cont.y=num.actives, max.per.page=num.actives)
        #         else
        #           hit.map(preds[[i]][,-1], y=y, y.labels=row.names(desc), descriptors=desc,
        #                   title=titles[i], max.cont.y=num.actives, max.per.page=num.actives, labels=names(preds[[i]])[-1])

      }
    }

    # Plot all descriptors for each method
    meths <- c()
    for (i in length(titles)) meths <- c(meths, names(preds[[i]])[-1])
    meths <- unique(meths)
    if (classify == "Y") {
      pmeths <- c()
      for (i in length(titles)) pmeths <- c(pmeths, names(probs[[i]])[-1])
      pmeths <- unique(pmeths)
      meths <- meths[!(meths %in% pmeths)]
      for (j in 1:length(pmeths)) {
        p <- data.frame()
        p.labels <- c()
        for (i in 1:length(titles)) {
          if (sum(names(probs[[i]]) %in% pmeths[j]) == 1) {
            if (ncol(p) == 0)
              p <- data.frame(probs[[i]][, pmeths[j]])
            else p <- data.frame(p, probs[[i]][, pmeths[j]])
            names(p)[ncol(p)] <- abrev.titles[i]
            p.labels <- c(p.labels, abrev.titles[i])
          }
        }
        if (ncol(p)>0) {
          HitCurve(p, y=y, title=paste("Split", splidx, ":", pmeths[j]), phat.labels=p.labels)
#           if (ncol(desc)==0)
#             hit.map(p, y=y, title=pmeths[j], cexCol=1, labels=p.labels)
#           else
#             hit.map(p, y=y, y.labels=row.names(desc), descriptors=desc, title=pmeths[j], cexCol=1, labels=p.labels)
#         }
        }
      }
    }
    for (j in 1:length(meths)) {
      p <- data.frame()
      p.labels <- c()
      for (i in 1:length(titles)) {
        if (sum(names(preds[[i]]) %in% meths[j]) == 1) {
          if (ncol(p) == 0)
            p <- data.frame(preds[[i]][, meths[j]]) else p <- data.frame(p, preds[[i]][, meths[j]])
            names(p)[ncol(p)] <- abrev.titles[i]
            p.labels <- c(p.labels, abrev.titles[i])
        }
      }
      if (ncol(p) > 0) {
        ContCurve(p, y = y, title = paste("Split", splidx, ":", meths[j]),
                  yhat.labels = p.labels)
        #         if (ncol(desc)==0)
        #           hit.map(p, y=y, title=meths[j],
        #                   cexCol=1, max.cont.y=num.actives, max.per.page=num.actives, labels=p.labels)
        #         else
        #           hit.map(p, y=y, y.labels=row.names(desc), descriptors=desc, title=meths[j],
        #                   cexCol=1, max.cont.y=num.actives, max.per.page=num.actives, labels=p.labels)
      }
    }

    #     #Create plot.pdf
    #     i<-length(titles)
    #     if ((classify=="Y") && (ncol(probs[[i]])>1)) {
    #       hit.curve(probs[[i]][,-1], y=y, title=titles[i], phat.labels=names(probs[[i]])[-1])
    #       cont.curve(preds[[i]][,!(names(preds[[i]]) %in% names(probs[[i]]))],
    #                  y=y, curves.only=TRUE, start.col=(ncol(probs[[i]])-1), title=titles[i],
    #                  yhat.labels=names(preds[[i]])[!(names(preds[[i]]) %in% names(probs[[i]]))])
    #     }
    #     else {
    #       cont.curve(preds[[i]][,-1], y=y, title=titles[i], yhat.labels=names(preds[[i]])[-1])
    #     }
  }
  # dev.off()

  #     #Create combined predictions
  #     head<-","
  #     p<-data.frame()
  #     for (i in 1:length(titles)) {
  #       head<-paste(head,",",titles[i],sep="")
  #       if (ncol(preds[[i]])>2)
  #         for (j in 2:(ncol(preds[[i]])-1))
  #           head<-paste(head,",",sep="")
  #         if (ncol(p)==0) p<-data.frame(preds[[i]])
  #         else p<-data.frame(p,subset(preds[[i]],select=-1),check.names=FALSE)
  #     }
  #     write.table(head,file=paste(root,dest,"/pred.csv",sep=""),col.names=FALSE,quote=FALSE,row.names=FALSE)
  #     IDS <- rownames(p)
  #     p <- cbind(IDS,p)
  #     write.table(p,file=paste(root,dest,"/pred.csv",sep=""),quote=FALSE,row.names=F,sep=",",append=TRUE)
  #
  #     #Create combined probabilities
  #     if (classify=="Y") {
  #       head<-","
  #       p<-data.frame()
  #       for (i in 1:length(titles)) {
  #         head<-paste(head,",",titles[i],sep="")
  #         if (ncol(probs[[i]])>2)
  #           for (j in 2:(ncol(probs[[i]])-1))
  #             head<-paste(head,",",sep="")
  #           if (ncol(p)==0) p<-data.frame(probs[[i]])
  #           else p<-data.frame(p,subset(probs[[i]],select=-1),check.names=FALSE)
  #       }
  #       write.table(head,file=paste(root,dest,"/prob.csv",sep=""),col.names=FALSE,quote=FALSE,row.names=FALSE)
  #       IDS <- rownames(p)
  #       p <- cbind(IDS,p)
  #       write.table(p,file=paste(root,dest,"/prob.csv",sep=""),quote=FALSE,row.names=FALSE,sep=",",append=TRUE)
  #     }


  #     #Create combined summaries
  #     sums<-c()
  #     if (bndir!="") sums<-c(sums,readLines(paste(root,bndir,"/summary.txt",sep="")))
  #     if (phdir!="") sums<-c(sums,readLines(paste(root,phdir,"/summary.txt",sep="")))
  #     if (apdir!="") sums<-c(sums,readLines(paste(root,apdir,"/summary.txt",sep="")))
  #     if (fpdir!="") sums<-c(sums,readLines(paste(root,fpdir,"/summary.txt",sep="")))
  #     if (chdir!="") sums<-c(sums,readLines(paste(root,chdir,"/summary.txt",sep="")))
  #     writeLines(sums,paste(root,dest,"/summary.txt",sep=""))
  # }
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
  # bitmap(paste(plotfile,'.bmp',sep=''))
  # pdf(paste(plotfile,'.pdf',sep=''))
  if ((classify == "Y") && (ncol(probs) > 1)) {
    HitCurve(probs[, -1], y = probs[, 1], phat.labels = names(probs)[-1])
    ContCurve(preds[, !(names(preds) %in% names(probs))],
              yhat.labels = names(preds)[!(names(preds) %in%
                                             names(probs))], y = preds[, 1],
              curves.only = TRUE, start.col = (ncol(probs) - 1))
  } else {
    ContCurve(preds[, -1], y = preds[, 1], yhat.labels = names(preds)[-1])
  }
  # dev.off()
  # pdf(paste(plotfile,'s.pdf',sep=''))
  if ((classify == "Y") && (ncol(probs) > 1)) {
    HitCurve(probs[, -1], y = probs[, 1], phat.labels = names(probs)[-1])
    ContCurve(preds[, !(names(preds) %in% names(probs))],
              yhat.labels = names(preds)[!(names(preds) %in%
                                             names(probs))], y = preds[, 1],
              curves.only = TRUE, start.col = (ncol(probs) - 1))
    # HitMap(data.frame(probs[,-1],preds[,!(names(preds) %in% names(probs))]),
    # y=preds[,1], y.labels=row.names(desc), descriptors=desc,
    # labels=c(names(probs)[-1],names(preds)[!(names(preds) %in% names(probs))]))
  } else {
    ContCurve(preds[, -1], y = preds[, 1], yhat.labels = names(preds)[-1])
    # HitMap(preds[,-1], y=preds[,1], y.labels=row.names(desc), descriptors=desc,
    # labels=names(preds)[-1])
  }
  # dev.off()
}
