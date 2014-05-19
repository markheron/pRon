\name{cut_out_fasta_multiple_from_one_chr}
\alias{cut_out_fasta_multiple_from_one_chr}
\title{cut_out_fasta_multiple_from_one_chr}
\usage{
cut_out_fasta_multiple_from_one_chr(pos, strand, size, order, chr_fasta)
}
\arguments{
  \item{pos}{vector of central position for cuting out the
  region}

  \item{strand}{vector of "+"/"-" to specify if the forward
  or the reverse complement sequence shall be used}

  \item{size}{region size before and after the positions to
  cut out}

  \item{order}{what oligonucleotide order will later be
  used (to extend the cut region for symetry)}

  \item{chr_fasta}{DNAString of the chromosome from which
  to cut out the regions}
}
\value{
DNAStringSet of cut out fragment
}
\description{
Cut's out multiple regions from a single DNAString (i.e.
typically chromosome).
}
\author{
Mark Heron
}

