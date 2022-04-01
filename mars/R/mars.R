# Lab 1, week 2
#Shell implementation of MARS algorithm
#student: Dawu Liu
#SFU ID: 301116278


#' Title
#'
#' @param formula a R formula
#' @param data a data
#' @param control a mars.control object
#'
#' @return
#' @export
#'
#' @examples
#' @import stats
mars <- function(formula, data, control=NULL, ...) {
  cc <- match.call()
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)
  x <- x[,-1]

  fwd <- fwd_stepwise(y, x, control)
  bwd <- bwd_stepwise(fwd, control)

  dat <- data.frame(y=bwd$y,bwd$B)
  model = lm(formula, dat)

  output <- list()
  output$call <- match.call()
  output$formula <- formula
  output$y <- bwd$y
  output$B <- bwd$B
  output$Bfuncs <- bwd$Bfuncs
  output$x_names <- colnames(x)
  output$model = model

  class(output) = "mars"
  return (output)
}

#' Title
#'
#' @param y parameter
#' @param x parameter
#' @param control parameter
#'
#' @return
#' @export
#'
#' @examples
fwd_stepwise <- function(y,x,control){
  if(control$Mmax<2) {
    warning("Mmax is not greater than or equal to 2")
    control$Mmax = 2
  }

  # Initialize:
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  B <- init_B(N,control$Mmax) # Exercise: write init_B()
  Bfuncs <- vector(mode = "list", length = control$Mmax+1)
  #---------------------------------------------------
  # Looping for forward selection:
  for(i in 1:(control$Mmax/2)) {
    M <- 2*i-1
    lof_best <- Inf
    for(m in 1:M) { # choose a basis function to split
      vv <- setdiff(1:n, Bfuncs[[m]][,2])#  variables v not already in Bfuncs[[m]], column 2 is v for the variable
      for(v in vv){ # select a variable to split on
        tt <- split_points(x[,v],B[,m]) # split points
        for(t in tt) {
          Bnew <- data.frame(B[,(1:M)],
                             Btem1=B[,m]*h(+1, x[,v], t),
                             Btem2=B[,m]*h(-1, x[,v], t))
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat,control) # updated LOF
          if(lof < lof_best) {
            lof_best <- lof
            split_best <- c(m=m,v=v,t=t)
          } # end if
        } # end loop over splits
      } # end loop over variables
    } # end loop over basis functions to split
    m <- split_best["m"]; v <- split_best["v"]; t <- split_best["t"]

    B[,M+1] <- B[,m]*h(-1, x[,v], t)
    B[,M+2] <- B[,m]*h(+1, x[,v], t)

    Bfuncs[[M+1]] <- rbind(Bfuncs[[m]], c(s=-1, v, t))
    Bfuncs[[M+2]] <- rbind(Bfuncs[[m]], c(s=+1, v, t))
  } # end loop over M
  colnames(B) <- paste0("B",(0:(ncol(B)-1)))
  return(list(y=y,B=B,Bfuncs=Bfuncs))
}



# take output of fwd_stepwise()
#' Title
#'
#' @param fwd parameter
#' @param control parameter
#'
#' @return
#' @export
#'
#' @examples
bwd_stepwise <- function(fwd, control) {
  Mmax <- control$Mmax
  Jstar <- 2:(Mmax+1)
  dat1 <- data.frame(y=fwd$y,fwd$B)
  lofstar <- LOF(y~.-1, dat1, control)
  Kstar <- Jstar
  for(M in (Mmax+1):2){
    L <- Kstar
    b <- Inf
    for(m in L){
      K <- setdiff(L, m)
      dat2 <- data.frame(y=fwd$y,fwd$B[,K])
      lof <- LOF(y~., dat2, control)
      if(lof < b){
        b <- lof
        Kstar <- K
      }
      if(lof < lofstar){
        lofstar <- lof
        Jstar <- K
      }
    }
  }
  Jstar <- c(1, Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
}

H <- function(x) {
  return(as.numeric(x>=0))
}

# pmax Returns the  parallel maxima of two vectors
h <- function(s,x,t){
  return(pmax(0, s*(x-t)))
}

LOF <- function(formula, data, control){
  N <- nrow(data)
  ff <- lm(formula, data)
  M <- ncol(data)-2 # number of non-constant basis functions
  Ctilde <- sum(hatvalues(ff)) + (control$d*M) # ˜C(M) = C(M)+dM
  RSS <- sum(residuals(ff)^2) # residual sum of squares
  GCV <- RSS * N/(N-Ctilde)^2
  return(GCV)
}

init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

# constructor
new_mars.control <- function(x = list()) {
  structure(x, class = "mars.control")
}

# validator
validate_mars.control <- function (x) {
  if(is.integer(x$Mmax) ==FALSE) warning("Mmax not integer")
  if(x$Mmax%%2 != 0) warning("Mmax not even")
  if(x$Mmax < 2) warning("Mmax is less than 2")
  if(is.numeric(x$d) == FALSE) warning("d is not numeric")
  if(is.logical(x$trace) == FALSE) warning("trace is not logical")
  return (out=list(Mmax=x$Mmax, d=x$d, trace=x$trace))
}

# helper
mars.control <- function(Mmax=2, d=3, trace=FALSE) {
  Mmax = as.integer(Mmax)
  x <- list(Mmax=Mmax, d=d, trace=trace)
  if (Mmax <2 ) {Mmax = 2}
  x <- validate_mars.control(x)
  new_mars.control(x)
}
