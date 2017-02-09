# removed option to provide y response vector, this will be used later in a
# different S3 ModelTrain method

ReadInData <- function(data, ycol, xcols, idcol, type = "ANALYZE") {
  
  if (is.na(idcol)) {
    if (!all(c(ycol, xcols) <=  ncol(data)))
      stop("there is a column number larger than number of columns in 'data'")
  } else {
    if (!all(c(ycol, xcols, idcol) <= ncol(data)))
      stop("there is a column number larger than number of columns in 'data'")
  }

  # might still need this for prediction

  #-----Read in data
  #   if (is.na(ycol)) missingy = TRUE
  #   else missingy = FALSE

  work.data.response <- subset(data, select = ycol)
  work.data <- subset(data, select = xcols)
  work.data <- cbind(work.data.response, work.data)
  rm(work.data.response)

  if (!is.na(idcol)) {
    row.names(work.data) <- data[, idcol]
  }
  names(work.data)[1] <- "y"

  rm.y <- is.na(work.data$y)
  if (sum(rm.y) > 0) {
    warning(paste("WARNING.....", sum(rm.y),
                  "missing responses were provided and will not be used"))
    work.data <- subset(work.data, !rm.y)
  }
  rm(rm.y)

  rm.desc <- apply(work.data, 1, function(x) sum(is.na(x)) > 0)
  if (sum(rm.desc) > 0) {
    warning(paste("WARNING.....", sum(rm.desc),
                  "rows with missing responses or descriptors were provided and will not be used"))
    work.data <- subset(work.data, !rm.desc)
  }
  rm(rm.desc)
  if (type == "ANALYZE") {
    if (length(unique(work.data$y)) == 1) {
      stop(paste("ERROR..... at least 2 different responses are required for analysis"))
    }
    # removing constant columns
    rm.cols <- (apply(work.data, 2, var) == 0)
    if (sum(rm.cols) > 0) {
      warning(paste("WARNING.....", sum(rm.cols),
                    "constant descriptor columns were provided and will not be used"))
      work.data <- subset(work.data, select = (!rm.cols))
    }
  } else rm.cols <- 0

  if (sum(!(apply(work.data, 1, is.numeric))) > 0) {
    stop(paste(
      "ERROR..... non-numeric responses or descriptors were found,
      responses and descriptors must be numeric",
      dim(work.data)[1]))
  }

  #   if (missingy) {
  #     work.data <- cbind( rep(NA,nrow(work.data)), work.data )
  #     names( work.data )[1] <- 'y'
  #   }

  return(list(work.data, rm.cols))
}

