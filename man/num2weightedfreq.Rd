\name{num2weightedfreq}
\alias{num2weightedfreq}
\title{num2weightedfreq}
\usage{
num2weightedfreq(seqnum, weights, oligo_length)
}
\arguments{
  \item{seqnum}{(numeric matrix) oligo-nucleotide-coded
  representation of fasta sequences}

  \item{weights}{of the diferent fasta sequences}

  \item{oligo_length}{(integer) length of the
  oligonucleotide-coding used (order+1)}
}
\value{
matrix of weighted oligonucleotide frequencies
}
\description{
Calculates the weighted mere frequencies given a
oligonucleotide-coded representation matrix
}

