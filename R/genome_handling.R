
#' get_genome_avg
#' 
#' Computes the genome-wide oligomer frequencies.
#' 
#' @export
#' @param genome_dir the directory where the chromosome fasta files can be found
#' @param oligo_length (numeric) length of the oligo-nucleotides to cound the average frequencies for
#' @return the frequencies of the different oligomers
#'
get_genome_avg <- function(genome_dir, oligo_length) {
  
  oligo_counts <- rep(0, 4^oligo_length)
  
  for(chr in maRs::get_files(genome_dir)) {
    
    sequence_as_num <- fasta2num( Biostrings::readDNAStringSet(paste(genome_dir,chr,sep="")), oligo_length)
    chr_counts <- table(sequence_as_num[sequence_as_num > 0])
    oligo_counts <- oligo_counts + chr_counts
  }
  return(oligo_counts/sum(oligo_counts))
}


#' read_genome_fasta
#' 
#' Load chromosome fasta files from a genome folder into a DNAStringSet object.
#' 
#' @export
#' @param genome_dir the directory where the chromosome fasta files can be found
#' @return DNAStringSet of the genome
#' 
read_genome_fasta <- function(genome_dir) {
  
  genome <- list()
  
  for(chr in maRs::get_files(genome_dir)) {
    genome[[sub("\\..*", "", chr)]] <- Biostrings::readDNAStringSet(paste(genome_dir,chr,sep=""))
  }
  genome <- Biostrings::DNAStringSet( lapply(genome, function (x) x[[1]]))
  return(genome)
}




##' cut_out_fasta_multiple
##'
##' Cuts out multiple regions from a genome (DNAStringSet).
##' 
##' @export
##' @param chr vector of chromosomes (DNAString) names where the sequences should be cutout from
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order will later be used (to extend the cut region for symetry)
##' @param genome_dir the directory where the chromosome fasta files can be found
##' @return DNAStringSet of cut out fragment
##' 
cut_out_fasta_multiple <- function(chr, pos, strand, size, order, genome_dir) {
  
  genome <- read_genome_fasta(genome_dir)
  
  seqs <- Biostrings::DNAStringSet(rep("",length(pos)))
  
  for(i in seq_along(seqs)) {
    
    seqs[[i]] <- cut_out_fasta_single(genome[[chr[i]]], start=pos[i]-size -((strand[i]=="-")*order), end=pos[i]+size +((strand[i]=="+")*order), strand[i] )
  }
  names(seqs) <- ""
  return(seqs)
}


##' cut_out_fasta_single
##' 
##' Cuts out single region from a single DNAString.
##' 
##' @export
##' @param fasta DNAString to cut out from
##' @param start of region to cut out
##' @param end of region to cut out
##' @param strand cut out forward or reverse complement
##' @return DNAString of cut out fragment
##' 
cut_out_fasta_single <- function(fasta, start, end, strand) {
  
  fasta_seq <- XVector::subseq(fasta, start, end)
  if(strand == "-") {
    fasta_seq <- Biostrings::reverseComplement(fasta_seq)
  }
  return(fasta_seq)
}


##' cut_out_fasta_multiple_from_one_chr
##'
##' Cuts out multiple regions from a single DNAString(Set) (i.e. typically a chromosome).
##' 
##' @export
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order will later be used (to extend the cut region for symetry)
##' @param chr_fasta DNAString(Set) of the chromosome from which to cut out the regions
##' @return DNAStringSet of cut out fragments
##' 
cut_out_fasta_multiple_from_one_chr <- function(pos, strand, size, order, chr_fasta) {
  
  seqs <- Biostrings::DNAStringSet(rep("",length(pos)))
  
  for(i in 1:length(seqs)) {
    
    seqs <- XVector::subseq(rep(Biostrings::DNAStringSet(chr_fasta), length(pos)), start=pos-size -((strand=="-")*order), end=pos+size +((strand=="+")*order))
    
  }
  names(seqs) <- ""
  return(seqs)
}



##' convertSparse2Complete_ff
##'
##' Converts sparse representation to complete ff vector representation for a list of data.
##' 
##' @export
##' @param sparse list of sparse representation matricies
##' @param lengths list of vector lengths, names must match those of \code{sparse}
##' @return list of ff vectors, which are the complete representation of each of \code{sparse} elements
##' 
convertSparse2Complete_ff <- function(sparse, lengths) {
  
  complete <- list()
  
  for(chr in names(sparse)) {
    complete[[chr]] <- ff::ff(0, length=lengths[[chr]])
    complete[[chr]][sparse[[chr]][,1]] <- sparse[[chr]][,2]
  }
  invisible(complete)
}




##' getRunningWindowCG
##'
##' Computes a running window C+G\% for a DNAStringSet.
##' 
##' @export
##' @param fastas (DNAStringSet) for which to calculate the running window C+G\%
##' @param half_window_size how much to extend the window in either direction
##' @return list of ff vectors that contain the running C+G\%
##' 
getRunningWindowCG <- function(fastas, half_window_size=73) {
  
  genome_list <- list()
  for(chr in names(fastas)) {
    seqnum <- fasta2num(fastas[chr],1)
    cg <- as.numeric((seqnum == 2) | (seqnum == 3))
    genome_list[[chr]] <-  ff::as.ff(maRs::smear(cg, from=-half_window_size, to=half_window_size)/(2*half_window_size+1))
  }
  invisible(genome_list)
}
