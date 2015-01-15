

##' fasta2sparse_boolean
##'
##' Creates sparse coding matrix of a fasta sequence.
##' @export
##' @param fasta DNAString
##' @param oligo_length oligonucleotide length to code the sequence in
##' @return sparse representation matrix of the DNAString
##' @author Mark Heron
fasta2sparse_boolean <- function(fasta, oligo_length) {
  
  fasta_num <- fasta2num(fasta, oligo_length)
  
  fasta_sparse <- matrix(FALSE, nrow=length(oligo_names(oligo_length)) ,ncol=length(fasta_num))
  
  for(i in 1:nrow(fasta_sparse)) {
    fasta_sparse[i,] <- fasta_num==i
  }
  #fasta_sparse[is.na(fasta_sparse)] <- FALSE
  
  return( fasta_sparse )
}




##' fasta2sparse_boolean_ff
##'
##' Creates sparse coding ff object of a fasta sequence.
##' @export
##' @param fasta DNAString
##' @param oligo_length oligonucleotide length to code the sequence in
##' @return sparse ff.matrix representation of the DNAString
##' @author Mark Heron
fasta2sparse_boolean_ff <- function(fasta, oligo_length) {
  
  fasta_num <- fasta2num(fasta, oligo_length)
  
  fasta_sparse <- ff(FALSE, vmode="boolean", dim=c(length(oligo_names(oligo_length)), length(fasta_num)))
  
  for(i in 1:nrow(fasta_sparse)) {
    fasta_sparse[i,] <- fasta_num==i
  }
  
  return( fasta_sparse )
}
