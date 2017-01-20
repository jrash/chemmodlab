
PrintTime <- function(pt, st) {
  cat("Real time used: ", (proc.time() - pt)[3], "\n")
  cat(" CPU time used: ", st[3], "\n\n\n")
}

UpdateStatus <- function(statusfile, text) {
  if ((!(is.na(statusfile))) && (nchar(statusfile) > 0)) {
    writeLines(text, statusfile)
  }
}
