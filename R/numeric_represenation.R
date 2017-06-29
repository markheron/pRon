

#' fasta2num
#' 
#' Transforms a DNAStringSet with equal length objects into a matrix of numeric representation of the oligo-nucleotides.
#' Any oligo-nucleotide with a letter differing from A,C,G,T is represented as **NA**
#' 
#' @export
#' @param fastas (DNAStringSet) of same length fasta sequences that should be converted (single DNAString also possible)
#' @param oligo_length (numeric) length of the oligo-nucleotides to use.
#' @param method ("lookupTable") one of"as.integer", "lookupTable", "slidingView" and "memoryLimited" they are different implementations that have different speed/memory tradeoffs
#' @param simplify (bool) in the case of a single fasta sequence, should a vector be returned instead of a matrix
#' @return matrix (or vector) of numeric representation of the fastas
#' 
fasta2num <- function (fastas, oligo_length, method="lookupTable", simplify=TRUE) {
  
  if(class(fastas) == "DNAString") {
    fastas <- Biostrings::DNAStringSet(fastas)
  }
  
  # could use chartr instead of this, but I have to strsplit them anyway, so this probably has similar speed
  base2num <- c("A"=1, "C"=2, "G"=3, "T"=4, "N"=NA)
  
  xinteger2num <- c(1,2,NA,3,NA,NA,NA,4)
  if(method == "as.integer") { # kind of a hack, numbering could change in the future
    seqnum_zero <- sapply(fastas, function(x) xinteger2num[as.integer(x)])
    
  } else if(method == "slidingView") {
    seqnum_zero <- sapply(fastas, function (x) colSums(t(letterFrequencyInSlidingView(x, 1, c("A","C","G","T")))*(1:4)))
    
  } else if (method == "lookupTable") { # previous line implementation is faster, see benchmark  (on human sized fasta's this seems faster and has less memory footprint)
    seqnum_zero <- sapply(fastas, function (fasta) base2num[strsplit(as.character(fasta),"")[[1]]])
    # seqnum_zero <- sapply(strsplit(as.character(fastas),""), function (fasta) base2num[fasta]) # similar speed
    
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
  rownames(seqnum_zero) <- NULL # because they don't make sense

  seqnum_zero[seqnum_zero == 0] <- NA
  
  seqnum_higher <- seqnum_zero
  if(oligo_length > 1) {
    for(i in 2:oligo_length) {
      seqnum_higher <- as.matrix( (seqnum_higher[-nrow(seqnum_higher),]-1)*4 + seqnum_zero[i:nrow(seqnum_zero),] )
    }
  }
  if(simplify & ncol(seqnum_higher)==1) {
    seqnum_higher <- as.vector(seqnum_higher)
  }
  
  return(seqnum_higher)
}

  



#' evenfasta2num
#' 
#' Transforms a DNAStringSet into a list of numeric representation of the mono-nucleotides.
#' Any mono-nucleotide that isn't A,C,G or T is represented as **NA**.
#' This is a streamlined version of \code{\link{fasta2num}} for mono-nucleotides that additionally allows varying sequence lengths.
#'
#' @export
#' @param fastas (DNAStringSet) that can be different length DNA sequences that should be converted (single DNAString also possible)
#' @return list of numeric representation of the DNA sequences
#'
evenfasta2num <- function (fastas) {
  
  # could use chartr instead of this, but I have to strsplit them anyway, so this probably has similar speed
  base2num <- c("A"=1, "C"=2, "G"=3, "T"=4, "N"=NA)
  
  seqnum_zero <- lapply(strsplit(as.character(fastas),""), function (x) base2num[x])
  
  return(seqnum_zero)
}





#' num2freq
#'
#' Calculates the oligo-nucleotide frequencies given a numerically encoded matrix representation.
#' 
#' @export
#' @param seqnum (numeric matrix) numerically encoded oligo-nucleotide representation of DNA sequences
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
#' Calculates the weighted oligo-nucleotide frequencies given a numerically encoded matrix representation.
#' 
#' @export
#' @inheritParams num2freq
#' @param weights (numeric) to be given to the fasta sequences (should have length matching \code{ncol(seqnum)})
#' @return matrix of weighted oligonucleotide frequencies
#' 
num2weightedfreq <- function(seqnum, weights, oligo_length) {
  
  oligos <- oligo_names(oligo_length)
  freqs <- sapply(1:length(oligos), function (oligo) colSums(apply((seqnum==oligo), 1, function (x) x*weights), na.rm=TRUE))
  colnames(freqs) <- oligos
  return(t(freqs))
}
