

#' fasta2num
#' 
#' Transforms a DNAStringSet with equal length objects into a matrix of numeric representation of the oligo-nucleotides.
#' Any oligo-nucleotide with a letter differing from A,C,G,T is represented as **NA**
#' 
#' @export
#' @param fastas (DNAStringSet) of same length fasta sequences that should be converted
#' @param oligo_length (numeric) length of the oligo-nucleotides to use.
#' @return matrix of numeric representation of the fastas
#' 
fasta2num <- function (fastas, oligo_length, method="lookupTable") {
  
  # could use chartr instead of this, but I have to strsplit them anyway, so this probably has similar speed
  base2num <- c("A"=1, "C"=2, "G"=3, "T"=4, "N"=NA)
  
  if(method == "slidingView") {
    seqnum_zero <- sapply(fastas, function (x) colSums(t(letterFrequencyInSlidingView(x, 1, c("A","C","G","T")))*(1:4)))
  } else if (method == "lookupTable") { # previous line implementation is faster, see benchmark  (on human sized fasta's this seems faster and has less memory footprint)
    seqnum_zero <- sapply(fastas, function (fasta) base2num[strsplit(as.character(fasta),"")[[1]]])
  } else if (method == "memoryLimited") {
      
    seqnum_zero <- sapply(fastas, function (fasta) {
      positions <- c( seq(1,length(fasta),by=10000000),length(fasta)+1)
      tmp <- integer(length(fasta))
      for( i in 1:(length(positions)-1)) {
        tmp[positions[i]:(positions[i+1]-1)] <- base2num[strsplit(as.character(substr(fasta,positions[i], positions[i+1]-1)),"")[[1]]]
      }
      return( tmp )
    })
  }

  seqnum_zero[seqnum_zero == 0] <- NA
  
  seqnum_higher <- seqnum_zero
  if(oligo_length > 1) {
    for(i in 2:oligo_length) {
      seqnum_higher <- as.matrix( (seqnum_higher[-nrow(seqnum_higher),]-1)*4 + seqnum_zero[i:nrow(seqnum_zero),] )
    }
  }
  return(seqnum_higher)
}




#' num2freq
#'
#' Calculates the mere frequencies given a oligo-nucleotide-coded representation matrix
#' @export
#' @param seqnum (numeric matrix) oligo-nucleotide-coded representation of fasta sequences
#' @param oligo_length (integer) length of the oligonucleotide-coding used (order+1)
#' @return matrix of oligonucleotide frequencies
#' 
num2freq <- function(seqnum, oligo_length) {
  
  oligos <- oligo_names(oligo_length)
  freqs <- sapply(1:length(oligos), function (oligo) rowMeans(seqnum==oligo, na.rm=TRUE))
  colnames(freqs) <- oligos
  return(t(freqs))
}


#' num2weightedfreq
#'
#' Calculates the weighted mere frequencies given a oligonucleotide-coded representation matrix
#' @export
#' @param seqnum (numeric matrix) oligo-nucleotide-coded representation of fasta sequences
#' @param weights of the different fasta sequences
#' @param oligo_length (integer) length of the oligonucleotide-coding used (order+1)
#' @return matrix of weighted oligonucleotide frequencies
#' 
num2weightedfreq <- function(seqnum, weights, oligo_length) {
  
  oligos <- oligo_names(oligo_length)
  freqs <- sapply(1:length(oligos), function (oligo) colSums(apply((seqnum==oligo), 1, function (x) x*weights), na.rm=TRUE))
  colnames(freqs) <- oligos
  return(t(freqs))
}




#' cut_out_seqnums_from_one_chr
#'
#' Cut's out multiple regions from a single oligonucleotide-coded respresentation vector (i.e. typically chromosome).
#' export not for now
#' @param pos vector of central position for cuting out the region
#' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
#' @param size region size before and after the positions to cut out
#' @param order what oligonucleotide order was used to code the fasta sequence
#' @param chr_num oligonucleotide-coded respresentation of the chromosome from which to cut out the regions
#' @return matrix of cut out fragment
#' 
cut_out_seqnums_from_one_chr <- function(pos, strand, size, order, chr_num) {
  
  # does this really do what it should do, if yes, why seqs should be over written and it isn't devined to start out with...?
  for(i in 1:length(seqs)) {
    
    seqs <-  vapply(pos, function (i) chr_num[i+(-size -((strand=="-")*order)):(size +((strand=="+")*order))], FUN.VALUE=rep(0, length((-size -((strand=="-")*order)):(size +((strand=="+")*order)))))
  }
  return(seqs)
}
