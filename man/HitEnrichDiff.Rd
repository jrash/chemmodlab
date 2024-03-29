% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hit_enrich_curve.R
\name{HitEnrichDiff}
\alias{HitEnrichDiff}
\title{Plot differences between hit enrichment curves}
\usage{
HitEnrichDiff(
  S.df,
  labels = NULL,
  y,
  x.max = NULL,
  log = TRUE,
  title = "",
  conf.level = 0.95,
  method = "sup-t",
  plus = TRUE,
  band.frac = NULL,
  yrange = NULL
)
}
\arguments{
\item{S.df}{Data frame where variables are numeric scores from at least 2
different algorithms. Rows represent unique compounds.}

\item{labels}{Character vector of labels for the different algorithms in
\code{S.df}. If missing, variable names in \code{S.df} will be used.}

\item{y}{Numeric vector of activity values. Activity values must be either 0
(inactive/undesirable) or 1 (active/desirable); no other values are
accepted. Compounds are assumed to be in the same order as in \code{S.df}.}

\item{x.max}{Integer, the maximum number of tests allowed on the x axis.}

\item{log}{Logical. \code{TRUE} plots the x axis on a log scale.}

\item{title}{Character string}

\item{conf.level}{Numeric, confidence coefficient}

\item{method}{Character indicates the method used to obtain confidence bands.
The default is \code{sup-t} but other options (not recommended) are
"theta-proj" and "bonf".}

\item{plus}{Logical. \code{TRUE} uses plus-adjusted version of \code{method}.}

\item{band.frac}{Numeric vector of fractions tested to be used in obtaining
confidence bands. Vector should be no longer than \code{y}, and should have
at least 20 entries. Entries should be in (0,1]. It is recommended that
entries be consistent with between 1 and \code{x.max} tests.}

\item{yrange}{Numeric vector of length 2. The desired range for the y axis.}
}
\description{
Plot differences between hit enrichment curves based on scores from multiple
algorithms. Actual activities are required. Additionally plot simultaneous
confidence bands for these differences. Plots may be used to determine if one
algorithm is "better" than another algorithm.
}
\details{
By default, \code{x.max} is \code{length(y)}, so that hit enrichment
curves are obtained for all observable fractions, i.e., fractions of
\code{(1:length(y))/length(y)}. By default, confidence bands are evaluated
based on a smaller grid of 40 fractions. This smaller grid is evenly spaced
on either the original grid of \code{(1:length(y))/length(y)}, or the log
scale of the original grid.
}
