#' Multivariate Adaptive Regression Splines (MARS)
#'
#' @description MARS (Multivariate Adaptive Regression Splines) is a regression analysis that automatically splits and better linear fit into non-linear models.
#'
#' @usage mars(formula, data, control=NULL, ...)
#'
#' @param formula an R formula
#' @param data a data frame containing the data
#' @param control an object of class 'mars.control'
#' @param ... further arguments
#' @details
#'  MARS - forward stepwise algorithm:
#'  Forward stepwise fits the data nicely by adjusting the coefficient values as well as deriving a proper set of basis functions. Forward stepwise produces basis functions which do not have zero pairwise product expectations. The advantage of this is that basis functions, except for B0, can be removed without leaving a hole in the predictor space.
#'
#'  MARS - backward stepwise algorithm:
#' Jstar(J*) backward stepwise is all of the basis function set and derived from forward stepwise. One basis function is deleted at each iteration of the outer loop, while the inner loop chooses which one to be deleted. The chosen basis function to be deleted is either the one improving the model by its removal or degrading the model the least by the removal. Constant basis function B0(x) = 1 is not considered for elimination. A sequence of (Mmax - 1) models is what is constructed from backward stepwise, and the current sequence has one less basis function than that of the previous sequence. Once the iteration is terminated, the users are left with the best model.
#'
#' @return an S3 object of class mars that includes the final regression model and a description of the basis functions which are constructed by the hinge functions
#'
#' @family methods for mars object
#' @export
#' @author Siyul Sam Byun, So Yeon Park, Dawu Liu
#' @references Jerame H. Friedman. "Multivariate Adaptive Regression Splines."
#' Ann, Statist. 19 (1) 1 - 67, March, 1991 \url{https://doi.org/10.1214/aos/1176347963}
#'
#' @seealso [plot.mars] for plotting the basis functions of the mars object
#' @seealso [predict.mars] for making predictions on a new data using the mars object
#' @seealso [summary.mars] for summarizing the mars object
#' @seealso [print.mars] for printing the mars object
#'
#' @examples test <- mars(y~x1+x2, data=mars::marstestdata, control=mars.control(Mmax=2))
#' @import stats
#'
mars <- function(formula, data, control=NULL, ...) {
  cc <- match.call()
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]
  if(is.null(control)) control <- mars.control()
  control <- validate_mars.control(control)

  fwd <- fwd_stepwise(y, x, control)
  bwd <- bwd_stepwise(fwd, control)

  model = lm(formula = y ~.-1, data = data.frame(y = y, bwd$B))

  output <-c(list(call=cc, formula=formula, y=y, B=bwd$B, Bfuncs=bwd$Bfuncs,
                 x_names=colnames(x)), model)

  class(output) = c("mars", class(model))
  return (output)
}

fwd_stepwise <- function(y,x,control){
  if(control$Mmax<2) {
    warning("Mmax is not greater than or equal to 2")
    control$Mmax = 2L
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

init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

bwd_stepwise <- function(fwd, control) {
  #Mmax <- control$Mmax
  Mmax <- ncol(fwd$B)-1
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

# pmax Returns the  parallel maxima of two vectors
h <- function(s,x,t){
  return(pmax(0, s*(x-t)))
}

LOF <- function(formula, data, control){
  N <- nrow(data)
  ff <- lm(formula, data)
  M <- length(coef(ff))-1 # number of non-constant basis functions
  Ctilde <- sum(diag(hatvalues(ff))) + (control$d*M) # ˜C(M) = C(M)+dM
  RSS <- sum(residuals(ff)^2) # residual sum of squares
  GCV <- RSS * N/(N-Ctilde)^2
  return(GCV)
}

split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}

# constructor, validator, and helper

new_mars.control <- function(control) {
  structure(control, class = "mars.control")
}

# validator
validate_mars.control <- function (control) {
  stopifnot(is.integer(control$Mmax), is.numeric(control$d), is.logical(control$trace))

  if(control$Mmax < 2) {
    warning("Mmax is less than 2, changing it to 2")
    control$Mmax <- 2L
  }

  if(control$Mmax%%2 != 0) {
    warning("Mmax not even, making it even")
    control$Mmax <- 2*ceiling(control$Mmax/2) # make Mmax even
    control$Mmax <- as.integer(control$Mmax)
  }

  return (control)
}


#' Constructor for 'mars.control' objects
#'
#' This function constructs a 'mars.control' object that specifies parameters used in the model fitting procedure.
#'
#' @param Mmax A maximum even integer of basis functions. Default value is 2.
#' @param d A smoothing parameter for the generalized cross-validation. Default value is 3.
#' @param trace A logical value
#'
#' @return a `mars.control` object
#' @export
#'
#' @examples test <- mars.control(Mmax = 10)
mars.control <- function(Mmax=2, d=3, trace=FALSE) {
  Mmax = as.integer(Mmax)
  x <- list(Mmax=Mmax, d=d, trace=trace)
  if (Mmax < 2) {Mmax = 2L}
  x <- validate_mars.control(x)
  new_mars.control(x)
}

