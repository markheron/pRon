## Functions to plot a DNA motif
## 
## Originally adapted by Holger from the seqLogo package of Bioconductor
## 
## name plot_logos
## author Holger Hartmann
## author Mark Heron
NULL


## get information content per position for pem
pwm2ic <- function(pwm) {
  
  ic <- 2 + colSums( pwm*log2(pwm), na.rm=TRUE) #na.rm=TRUE handles the case of 0 counts
  return(ic)
}

## get information gain per position for pem
pwm2ig <- function(pwm, bg) {
  
  rc <- colSums( pwm*log2(pwm/bg), na.rm=TRUE) #na.rm=TRUE handles the case of 0 counts
  return(rc)
}

######################
##
## plot sequence logo
##
######################

letterA <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  #x <- c(0, 4, 6,10, 8,6.8,3.2,2,0,  3.6,  5,6.4,3.6)
  #y <- c(0,10,10, 0, 0,3.0,3.0,0,0,  4.0,7.5,  4,  4)
  x <- c(0, 4, 6,10, 8,5.0,2,0,  3.2,3.6,6.4,6.8,3.2)
  y <- c(0,10,10, 0, 0,7.5,0,0,  3.0,4.0,4.0,3.0,3.0)
  
  x <- 0.1*x
  y <- 0.1*y
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
  }
  
  fill <- c("#008000","#008000")
  col <- c("#008000","#008000")
  
  list(x=x,y=y,id=id,fill=fill,col=col)
}

## T
letterT <- function(x.pos,y.pos,ht,wt,id=NULL){
  
  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- rep(1,8)
  }else{
    id <- rep(id,8)
  }
  
  fill <- "#FF0000"
  col <-  "#FF0000"
  
  list(x=x,y=y,id=id,fill=fill,col=col)
}

## C
letterC <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }
  
  fill <- "#0000FF"
  col <- "#0000FF"
  
  list(x=x,y=y,id=id,fill=fill,col=col)
}


## G
letterG <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)
  
  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)
  
  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))
  
  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)
  
  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)
  
  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)
  
  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))
  
  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  
  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)
  
  
  
  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
  }
  
  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)
  
  x <- x.pos + wt*x
  y <- y.pos + ht*y
  
  
  fill <- c("#FFA500","#FFA500")
  col  <- c("#FFA500","#FFA500")
  
  list(x=x,y=y,id=id,fill=fill,col=col)
  
}



addLetter <- function(letters,which,x.pos,y.pos,ht,wt){
  
  if (which == "A"){
    letter <- letterA(x.pos,y.pos,ht,wt)
  }else if (which == "C"){
    letter <- letterC(x.pos,y.pos,ht,wt)
  }else if (which == "G"){
    letter <- letterG(x.pos,y.pos,ht,wt)
  }else if (which == "T"){
    letter <- letterT(x.pos,y.pos,ht,wt)
  }else{
    stop("which must be one of A,C,G,T")
  }
  
  letters$x <- c(letters$x,letter$x)
  letters$y <- c(letters$y,letter$y)
  
  lastID <- ifelse(is.null(letters$id),0,max(letters$id))
  letters$id <- c(letters$id,lastID+letter$id)
  letters$fill <- c(letters$fill,letter$fill)
  letters$col <- c(letters$col,letter$col)
  letters
}

##' seqLogo
##' 
##' plot a sequence logo
##' 
##' @export
##' 
##' @param pwm position weight matrix to create a Sequence Logo from
##' @param ytype (character) the type of y-axis to be used "information" (or "ic", "content") for information content (default); "relative" (or "gain") for information gain; "probabilities" (or anything else) for probabilities
##' @param xaxis (logical) should an x-axis be drawn or vector of labels for the xaxis (must be same length as pwm)
##' @param yaxis (logical) should an y-axis be drawn
##' @param xfontsize (integer) size for x-axis text
##' @param yfontsize (integer) size for y-axis text
##' @param region in which to plot the figure, set by four margin values, set to \code{par("mar")} for default full plot
##' 
##' @examples
##' pwm <- rbind( A=c( 0,  0, 15,  5,  5,  1, 0.1),
##'               C=c(10,  0,  0,  5,  5,  2, 0.1),
##'               G=c(10,  0,  5,  0,  5,  3, 0.1),
##'               T=c( 0, 20,  0, 10,  5, 14,19.7)) /20
##'               
##' seqLogo(pwm, region=par("mar"))
##' 
##' @import grid
seqLogo <- function(pwm, background=c("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25), ytype="information", xaxis=TRUE, yaxis=TRUE, ylim=NULL, xlab="Position", xfontsize=15, yfontsize=15, region=c(0,0,1,1), ...){
  
  if (class(pwm) == "pwm"){
    pwm <- pwm@pwm
  }else if (class(pwm) == "data.frame"){
    pwm <- as.matrix(pwm)
  }else if (class(pwm) != "matrix"){
    stop("pwm must be of class pwm, matrix or data.frame")
  }
  
  if (any(abs(1 - colSums(pwm)) > 0.01)) {
    stop("Columns of PWM must add up to 1.0")
  }
  
  if (abs(1 - sum(background)) > 0.01) {
    stop("background must add up to 1.0")
  }
  
  chars <- c("A","C","G","T")
  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  
  if (ytype == "information" || ytype == "ic" || ytype == "content") {
    if(length(ylim) == 0) {
      ylim <- c(0,2)
    }
    ylab <- "Information content"
    facs <- pwm2ic(pwm)
  } else if (ytype == "relative" || ytype == "gain") {
    if(length(ylim) == 0) {
      ylim <- c(0,2)
    }
    ylab <- "Information gain"
    facs <- pwm2ig(pwm, background)
  } else {
    if(length(ylim) == 0) {
      ylim <- c(0,1)
    }
    ylab <- "Probability"
    facs <- rep(1,  ncol(pwm))
  }
  
  wt <- 1
  white_space <- 0.01
  
  x.pos <- 0
  for (j in 1: ncol(pwm)){
    
    hts <- pwm[,j]*facs[j]
    letterOrder <- order(hts)
    
    y.pos <- 0
    for (i in 1:4){
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht > 0) {
        if( (y.pos > 0) && (ht > white_space) ) {
          if(ht < 3*white_space) {
            y.pos <- y.pos + ht/3
            ht <- ht - ht/3
            letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt)
          } else {
            y.pos <- y.pos + white_space
            ht <- ht - white_space
            letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt)
          }
        } else {
          letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt)
        }
      }
      y.pos <- y.pos + ht
    }
    x.pos <- x.pos + wt
  }
  
  pushViewport(plotViewport(region, ...))
  pushViewport(dataViewport(xscale=c(0,ncol(pwm)),yscale=ylim,name="vp1", ...))
  
  grid.polygon(x=unit(c(0,ncol(pwm)+1,ncol(pwm)+1,0),"native"), y=unit(c(0,0,2,2),"native"),
               gp=gpar(fill="transparent",col="transparent"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id,
               gp=gpar(fill=letters$fill,col=letters$col,lwd=1))
  if(length(xaxis) == ncol(pwm) && length(xaxis) > 1) {
    grid.xaxis(at=seq(0.5,ncol(pwm)-0.5),label=xaxis, gp=gpar(fontsize=xfontsize))
    grid.text(xlab,y=unit(-3,"lines"), gp=gpar(fontsize=xfontsize))
  } else if (isTRUE(xaxis)) {
    grid.xaxis(at=seq(0.5,ncol(pwm)-0.5),label=1:ncol(pwm), gp=gpar(fontsize=xfontsize))
    grid.text(xlab,y=unit(-3,"lines"), gp=gpar(fontsize=xfontsize))
  }
  if (yaxis){
    grid.yaxis(gp=gpar(fontsize=yfontsize))
    grid.text(ylab,x=unit(-3,"lines"),rot=90, gp=gpar(fontsize=yfontsize))
  }
  popViewport()
  popViewport()
  par(ask=FALSE)
}

revComp <- function(pwm){
  reverse_pwm <- pwm # for correct col and row names
  reverse_pwm[] <- pwm[nrow(pwm):1, ncol(pwm):1]
  return(reverse_pwm)
}

