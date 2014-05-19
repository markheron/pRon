\name{cut_out_seqnums_from_one_chr}
\alias{cut_out_seqnums_from_one_chr}
\title{cut_out_seqnums_from_one_chr}
\usage{
cut_out_seqnums_from_one_chr(pos, strand, size, order, chr_num)
}
\arguments{
  \item{pos}{vector of central position for cuting out the
  region}

  \item{strand}{vector of "+"/"-" to specify if the forward
  or the reverse complement sequence shall be used}

  \item{size}{region size before and after the positions to
  cut out}

  \item{order}{what oligonucleotide order was used to code
  the fasta sequence}

  \item{chr_num}{oligonucleotide-coded respresentation of
  the chromosome from which to cut out the regions}
}
\value{
matrix of cut out fragment
}
\description{
Cut's out multiple regions from a single
oligonucleotide-coded respresentation vector (i.e.
typically chromosome). export not for now
}

