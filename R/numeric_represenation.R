

#' fasta2num
#' 
#' @export
#' @param fastas (DNAStringSet) of same length fasta sequences that should be converted
#' @param oligo_length (numeric) length of the oligo-nucleotides to use.
#' @return matrix of numeric representation of the fastas
#' 
fasta2num <- function (fastas, oligo_length) {
  
  # could use chartr instead of this, but I have to strsplit them anyway, so this probably has similar speed
  
  base2num <- c(1:4,NA)
  names(base2num) <- c("A","C","G","T", "N")
  
  seqnum_zero <- sapply(fastas, function (fasta) base2num[strsplit(as.character(fasta),"")[[1]]])
  
  seqnum_higher <- seqnum_zero
  if(oligo_length > 1) {
    for(i in 2:oligo_length) {
      seqnum_higher <- as.matrix((seqnum_higher[-dim(seqnum_higher)[1],]-1)*4+seqnum_zero[i:dim(seqnum_zero)[1],])
    }
  }
  return(seqnum_higher)
}




##' num2freq
##'
##' Calculates the mere frequencies given a oligonucleotide-coded representation matrix
##' @export
##' @param seqnum oligonucleotide-coded representation of fasta sequences
##' @param oligo_length length of the oligonucleotide-coding used (order+1)
##' @return matrix of oligonucleotide frequencies
##' @author Mark Heron
num2freq <- function(seqnum, oligo_length) {
  
  oligos <- oligo_names(oligo_length)
  freqs <- sapply(1:length(oligos), function (oligo) rowMeans(seqnum==oligo, na.rm=TRUE))
  colnames(freqs) <- oligos
  return(t(freqs))
}


##' num2weightedfreq
##'
##' Calculates the weighted mere frequencies given a oligonucleotide-coded representation matrix
##' @export
##' @param seqnum oligonucleotide-coded representation of fasta sequences
##' @param weights of the diferent fasta sequences
##' @param oligo_length length of the oligonucleotide-coding used (order+1)
##' @return matrix of weighted oligonucleotide frequencies
##' @author Mark Heron
num2weightedfreq <- function(seqnum, weights, oligo_length) {
  
  oligos <- oligo_names(oligo_length)
  freqs <- sapply(1:length(oligos), function (oligo) colSums(apply((seqnum==oligo), 1, function (x) x*weights), na.rm=TRUE))
  colnames(freqs) <- oligos
  return(t(freqs))
}




##' cut_out_seqnums_from_one_chr
##'
##' Cut's out multiple regions from a single oligonucleotide-coded respresentation vector (i.e. typically chromosome).
##' @export
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order was used to code the fasta sequence
##' @param chr_num oligonucleotide-coded respresentation of the chromosome from which to cut out the regions
##' @return matrix of cut out fragment
##' @author Mark Heron
cut_out_seqnums_from_one_chr <- function(pos, strand, size, order, chr_num) {
  
  for(i in 1:length(seqs)) {
    
    seqs <-  vapply(pos, function (i) chr_num[i+(-size -((strand=="-")*order)):(size +((strand=="+")*order))], FUN.VALUE=rep(0, length((-size -((strand=="-")*order)):(size +((strand=="+")*order)))))
    
  }
  return(seqs)
}
