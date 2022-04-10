#' Prints out the call and coefficients of a `mars` object
#'
#' @param x a `mars` object
#' @param ... further arguments
#'
#' @family methods
#' @export
#'
#' @examples mar <- mars(y~x1+x2+x3, data=mars::marstestdata, control=mars.control(Mmax=6))
#' print(mar)
#' @import stats
print.mars <- function (x, ...){
  cat("Call: \n ")
  print(x$call)
  cat("Coefficients: \n")
  print(x$coefficients)
  invisible(x)
}
