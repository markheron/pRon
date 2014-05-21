\name{plotOligoFreqs}
\alias{plotOligoFreqs}
\title{plotOligoFreqs}
\usage{
plotOligoFreqs(freqs, x_pos = NULL, ylim = NULL,
  xlab = "distance to position of interest", legend_nrow = 1, ...)
}
\arguments{
  \item{freqs}{(numeric matrix) oligonucleotide
  frequencies}

  \item{x_pos}{(numeric vector) x axis positions}

  \item{ylim}{(numeric tuple) range of the y axis}

  \item{xlab}{(character) x axis label, defaults to:
  \"distance to position of interest\"}

  \item{legend_nrow}{(integer) over how many rows should
  the legend be split}

  \item{...}{further parameters for \code{\link{plot}}}
}
\description{
plots oligonucleotide frequencies over a region as lines
}

