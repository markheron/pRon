% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_boolean_representation.R
\name{fasta2sparse_boolean_ff}
\alias{fasta2sparse_boolean_ff}
\title{fasta2sparse_boolean_ff}
\usage{
fasta2sparse_boolean_ff(fasta, oligo_length)
}
\arguments{
\item{fasta}{(DNAString) sequence to encode}

\item{oligo_length}{(integer) oligo-nucleotide length to encode the sequence in}
}
\value{
(boolean ff.matrix) sparse encoding of the DNAString, columns=positions, rows=oligo-nucleotides
}
\description{
Creates sparse encoding ff object of a DNA sequence.
}
