#' Predicted values based on mars object.
#'
#' @param object a `mars` object
#' @param newdata An optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted values are used.
#' @param ... further arguments
#'
#' @return predicted values of the response variable
#' @family methods
#' @export
#'
#' @examples mar <- mars(y~x1+x2+x3, data=dataset)
#'predict(object=mar, newdata=testdata)
predict.mars <- function(object, newdata, ...) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}

make_B <- function(X, Bfuncs){
  B <- data.frame(matrix(1,nrow=nrow(X),ncol=length(Bfuncs)))
  for (i in 2:(length(Bfuncs))){
    for (j in 1:nrow(Bfuncs[[i]])) {
      # B = h(s=, x[,v], t)*h(s=, x[,v], t)*h(s=, x[,v], t)..refer to lab 5 formula
      B[,i] <-  B[,i] * h(s=Bfuncs[[i]][j,][1], x=X[, Bfuncs[[i]][j,][2]] , t=Bfuncs[[i]][j,][3])
    }
  }
  B <- as.matrix(B)
  return (B)
}
