#' Summary method for class "mars".
#'
#' Prints function call, five-number summary, summary of hinge functions for each basis function,
#' and the coefficients of each basis function of a `mars` object
#'
#' @param object a `mars` object for which a summary is desired.
#' @param ... further arguments
#'
#' @family methods
#' @export
#'
#' @examples mar <- mars(y~x1+x2+x3, data=mars::marstestdata, control=mars.control(Mmax=6))
#' summary(mar)
#' @import stats
summary.mars <-function (object, ...){
  # Print Call and Residuals information
  five <- fivenum(object$residuals)
  names(five) <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("Call: \n")
  print(object$call)
  cat("Residuals: \n")
  print(five)
  cat("Residual degrees-of-freedom:", object$df.residual, "\n")

  # Basis function information
  bf <- object$Bfuncs
  names(bf) <- names(object$coefficients)
  cat("Hinge function: s(sign), v(variable), t(knot) \n")
  print(bf)

  cat("Coefficients: \n")
  print(object$coefficients)
}
