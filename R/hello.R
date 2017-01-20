# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!?")
}

library("htmltools")
writeNewBlogEntry <- function(curPkg, curVer, reposurl) {
  blogpost <-
    file.path("C:/Users/Vestige/Dropbox/chemmodlab_package/chemmodlab/utilities/description.txt")
  con <- file(blogpost, "wt")
  cat("New package", curPkg, "with initial version", curVer,"\n\n", file=con)
  dcf <- read.dcf("C:/Users/Vestige/Dropbox/chemmodlab_package/chemmodlab/DESCRIPTION")
  for (i in 1:ncol(dcf)) {
    cat("<strong>", colnames(dcf)[i], "</strong>: ",
        htmlEscape(dcf[1,i]), "<br>\n", sep="", file=con)
  }
  # closeBlogPost(con, reposurl, curPkg)
}

