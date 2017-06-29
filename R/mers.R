
#' create oligo-nucleotide names
#'
#' Creates an array with all possible oligo-nucleotides of the length \code{mer_length}, sorted alphabetically.
#' 
#' @export
#' @param oligo_length (numeric) length of the oligo-nucleotides to use.
#' @return (character vector) oligo-nucleotide name vector
#' @examples
#' oligo_names(2)
#' 
oligo_names <- function(oligo_length) {
  
  oligo_nuc_names <- c("")
  for(not_used in 1:oligo_length) {
    oligo_nuc_names <- paste0(rep(oligo_nuc_names, each=4), c('A','C','G','T'))
  }
  return(oligo_nuc_names)
}


#' complementary_oligo_positions
#'
#' Creates a position reference to the reverse complementary oligo-nucleotides in an alphabetically sorting.
#'
#' @export
#' @param oligo_length (numeric) length of the oligo-nucleotides to use.
#' @return (numeric vector) position of the matching reverse complement oligo-nucleotides
#' @examples
#' complementary_oligo_positions(2)
#' 
complementary_oligo_positions <- function(oligo_length) {
  
  swap_vector <- 0
  for(o in 1:oligo_length) {
    swap_vector <- rep(swap_vector, each=4) + (3:0)*(4^(o-1))
  }
  swap_vector <- swap_vector+1
  return(swap_vector)
}
