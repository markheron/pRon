
#' get_genome_avg
#' 
#' @export
#' @param genome_dir the directory where the chromosome fasta files can be found
#' @param oligo_length (numeric) length of the oligo-nucleotides to cound the average frequencies for
#' @return the frequencies of the different oligomers
#'
get_genome_avg <- function(genome_dir, oligo_length) {
  
  oligo_counts <- rep(0, 4^oligo_length)
  
  for(chr in get_files(genome_dir)) {
    
    sequence_as_num <- fasta2num( readDNAStringSet(paste(genome_dir,chr,sep="")), oligo_length)
    chr_counts <- table(sequence_as_num[sequence_as_num > 0])
    oligo_counts <- oligo_counts + chr_counts
  }
  return(oligo_counts/sum(oligo_counts))
}


#' read_genome_fasta
#' 
#' @export
#' @param genome_dir the directory where the chromosome fasta files can be found
#' @return list of DNAStringSets
#' 
read_genome_fasta <- function(genome_dir) {
  
  genome <- list()
  
  for(chr in get_files(genome_dir)) {
    genome[sub("\\..*", "", chr)] <- readDNAStringSet(paste(genome_dir,chr,sep=""))
  }
  genome <- DNAStringSet( lapply(genome, function (x) x[[1]]))
  return(genome)
}




##' cut_out_fasta_multiple
##'
##' Cut's out multiple regions from a genome (DNAStringSet).
##' @export
##' @param chr vector of chromosomes (DNAString) names where the sequences should be cutout from
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order will later be used (to extend the cut region for symetry)
##' @param genome_dir the directory where the chromosome fasta files can be found
##' @return DNAStringSet of cut out fragment
##' @author Mark Heron
cut_out_fasta_multiple <- function(chr, pos, strand, size, order, genome_dir) {
  
  genome <- read_genome_fasta(genome_dir)
  
  seqs <- DNAStringSet(rep("",length(pos)))
  
  for(i in 1:length(seqs)) {
    
    seqs[i] <- cut_out_fasta_single(genome[[chr[i]]], start=pos[i]-size -((strand[i]=="-")*order), end=pos[i]+size +((strand[i]=="+")*order), strand[i] )
  }
  names(seqs) <- ""
  return(seqs)
}


##' cut_out_fasta_single
##' 
##' Cut's out single region from a single DNAString.
##' @export
##' @param fasta DNAString to cut out from
##' @param start of region to cut out
##' @param end of region to cut out
##' @param strand cut out forward or reverse complement
##' @return DNAString of cut out fragment
##' @author Mark Heron
cut_out_fasta_single <- function(fasta, start, end, strand) {
  
  fasta_seq <- subseq(fasta, start, end)
  if(strand == "-") {
    fasta_seq <- reverseComplement(fasta_seq)
  }
  return(fasta_seq)
}


##' cut_out_fasta_multiple_from_one_chr
##'
##' Cut's out multiple regions from a single DNAString (i.e. typically chromosome).
##' @export
##' @param pos vector of central position for cuting out the region
##' @param strand vector of "+"/"-" to specify if the forward or the reverse complement sequence shall be used
##' @param size region size before and after the positions to cut out
##' @param order what oligonucleotide order will later be used (to extend the cut region for symetry)
##' @param chr_fasta DNAString of the chromosome from which to cut out the regions
##' @return DNAStringSet of cut out fragment
##' @author Mark Heron
cut_out_fasta_multiple_from_one_chr <- function(pos, strand, size, order, chr_fasta) {
  
  seqs <- DNAStringSet(rep("",length(pos)))
  
  for(i in 1:length(seqs)) {
    
    seqs <- subseq(rep(chr_fasta, length(pos)), start=pos-size -((strand=="-")*order), end=pos+size +((strand=="+")*order))
    
  }
  names(seqs) <- ""
  return(seqs)
}



##' convertSparse2Complete_ff
##'
##' Converts sparse representation to complete representation as ff vectors of genomic tags.
##' @export
##' @param sparse list of sparse representation matricies
##' @param lengths list of chromosome lengths, names must match those of sparse
##' @return list of ff vectors with complete representation of genomic tags
##' @author Mark Heron
convertSparse2Complete_ff <- function(sparse, lengths) {
  
  complete <- list()
  
  for(chr in names(sparse)) {
    complete[[chr]] <- ff(0, length=lengths[[chr]])
    complete[[chr]][sparse[[chr]][,1]] <- sparse[[chr]][,2]
  }
  invisible(complete)
}




##' get_dyad_pos
##'
##' Extracts the dyad positions in the chosen way from a table with mapped fragment start and ends
##' @export
##' @param data_list list of mapped fragments for each chromosome
##' @param dyad_base based on what position should the dyad position be calculated
##' @param offset offset if the dyad position isn't calculated from the center
##' @return list of ff matricies with dyad positions and intensity
##' @author Mark Heron
get_dyad_pos <- function(data_list, dyad_base="center", offset=73) {
  
  dyad_pos <- list()
  
  to_ffdf_table <- function (x) {
    return( as.ffdf(as.data.frame(table(x))) )
  }
  
  for(chr_name in names(data_list)) {
    
    if(dyad_base == "center") {
      tmp <- to_ffdf_table( floor((data_list[[chr_name]][,1] + data_list[[chr_name]][,2])/2) )
    } else if(dyad_base == "start") {
      tmp <- to_ffdf_table( data_list[[chr_name]][,1]+offset )
    } else if(dyad_base == "end") {
      tmp <- to_ffdf_table( data_list[[chr_name]][,2]-offset )
    } else if(dyad_base == "dinucleosome") {
      tmp <- to_ffdf_table( c(data_list[[chr_name]][,1]+offset, data_list[[chr_name]][,2]-offset) )
    }
    
    dyad_pos[[chr_name]] <- as.ff(matrix(c(as.numeric(as.character(tmp[,1])), tmp[,2]), ncol=2))
  }
  return(dyad_pos)
}



##' adjust_X_chr
##'
##' Adjusts lower X chromosome counts due cells being male or a mixture of male/female.
##' @export
##' @param ff_list list of ff objects each representing data of one chromosome
##' @param X_chr name of the X chromosome in ff_list
##' @param Xfactor factor by which the counts should be adjusted (2 for male celllines, 4/3 for male/female mixtures)
##' @return adjusted ff_list
##' @author Mark Heron
adjust_X_chr <- function(ff_list, X_chr, Xfactor) {
  
  ff_list[[X_chr]][,2] <- ff_list[[X_chr]][,2]*Xfactor
  return(ff_list)
}
