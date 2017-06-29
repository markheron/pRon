

##' fasta2sparse_boolean
##'
##' Creates sparse encoded matrix of a DNA sequence.
##' 
##' @export
##' @param fasta (DNAString) sequence to encode
##' @param oligo_length (integer) oligo-nucleotide length to encode the sequence in
##' @return (boolean matrix) sparse encoding of the DNAString, columns=positions, rows=oligo-nucleotides
##' 
fasta2sparse_boolean <- function(fasta, oligo_length) {
  
  fasta_num <- fasta2num(fasta, oligo_length)
  
  fasta_sparse <- matrix(FALSE, nrow=length(oligo_names(oligo_length)) ,ncol=length(fasta_num))
  
  for(i in 1:nrow(fasta_sparse)) {
    fasta_sparse[i,] <- fasta_num==i
  }
  
  return( fasta_sparse )
}




##' fasta2sparse_boolean_ff
##'
##' Creates sparse encoding ff object of a DNA sequence.
##' 
##' @export
##' @inheritParams fasta2sparse_boolean 
##' @return (boolean ff.matrix) sparse encoding of the DNAString, columns=positions, rows=oligo-nucleotides
##' 
fasta2sparse_boolean_ff <- function(fasta, oligo_length) {
  
  fasta_num <- fasta2num(fasta, oligo_length)
  
  fasta_sparse <- ff(FALSE, vmode="boolean", dim=c(length(oligo_names(oligo_length)), length(fasta_num)))
  
  for(i in 1:nrow(fasta_sparse)) {
    fasta_sparse[i,] <- fasta_num==i
  }
  
  return( fasta_sparse )
}
