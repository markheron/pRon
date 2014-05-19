

# #' plotGenomicCutouts
# #' 
# #' 
# plotGenomicCutouts <- function(chr, pos, strand, size, order, genome_folder) {
#   
#   fasta <- cut_out_fasta_multiple(chr, pos, strand, size, order, genome_folder)
#   
#   seqnum <- fasta2num(fasta, oligo_length=order+1)
#   
#   freqs <- num2freq(seqnum, order+1)
#   
#   plotOligoFreqs(freqs, x_pos=-size:size)
# }



#' plotOligoFreqs
#' 
#' @export
#' 
plotOligoFreqs <- function(freqs, x_pos=NULL, main="", ylim=NULL, legend_nrow=1) {
  
  if(all(is.nan(freqs))) {
    epsilons[] <- 0
  }
  
  if(length(ylim) == 0) {
    ylim <- c(min(freqs, na.rm =T),max(freqs, na.rm =T))
  }
  
  if(length(x_pos) == 0) {
    x_pos <- 1:dim(epsilons)[2]
  }
  
  colors <- distinctive_colors(dim(freqs)[1])
  
  plot(x_pos ,freqs[1,],typ='l',col=colors[1],ylim=ylim,xlab="", ylab="freq",xaxt='n', main=main)
  axis(side=1, pretty(x_pos, 10), mgp=c(3.5,1.5,0))
  title(xlab="distance to position of interest", mgp=c(4.5,1.5,0))
  for(i in 2:dim(freqs)[1]) {
    lines(x_pos ,freqs[i,],typ='l',col=colors[i])
  }
  legend("top",legend=rownames(freqs),fill=colors, ncol=dim(freqs)[1]/legend_nrow, cex=min(par('cex'), 6/(dim(freqs)[1]^(0.25))^2)) #,horiz=TRUE)
}
