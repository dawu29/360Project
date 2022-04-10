#' Plot a mars object
#'
#' Plot basis functions that are constructed with a single explanatory variable
#' or two explanatory variables
#'
#' @param x a mars object
#' @param ... further arguments
#'
#' @export
#'
#' @examples mar <- mars(y~x1+x2+x3, data=mars::marstestdata, control=mars.control(Mmax=6))
#'plot(x=mar)
#' @import graphics stats
plot.mars <- function(x, ...){
  data <- eval(x$call$data)
  tt <- terms(x$formula,data=data)
  tt <- delete.response(tt)
  mf <- model.frame(tt,data)
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)[,-1]
  Bf <- x$Bfuncs
  singleB <- which(sapply(Bf, function(x) NROW(x)==1))
  doubleB <- which(sapply(Bf, function(x) NROW(x)==2))
  n <- ceiling(sqrt(length(singleB)+length(doubleB)))
  opar <- par(mfrow=c(n,n),mar=c(2,2,2,2)); on.exit(par(opar))
  for(i in singleB){
    vv <- Bf[[i]][1,"v"]; varname1 <- x$x_names[[vv]]
    xx <- seq(from=min(X[,vv]),to=max(X[,vv]),length=100)
    yy <- h(x=xx,s=Bf[[i]][1,"s"],t=Bf[[i]][1,"t"])
    plot(xx,yy,type="l",xlab=varname1,main=varname1,...)
  }
  for(i in doubleB){
    vv1 <- Bf[[i]][1,"v"]; varname1 <- x$x_names[[vv1]]
    vv2 <- Bf[[i]][2,"v"]; varname2 <- x$x_names[[vv2]]
    xx <- seq(from=min(X[,vv1]),to=max(X[,vv1]),length=100)
    yy <- seq(from=min(X[,vv2]),to=max(X[,vv2]),length=100)
    ff <- function(x,y) {
      h(x=x,s=Bf[[i]][1,"s"],t=Bf[[i]][1,"t"])*
        h(x=y,s=Bf[[i]][2,"s"],t=Bf[[i]][2,"t"])}
    zz <- outer(xx,yy,FUN=ff)
    persp(xx,yy,zz,xlab=varname1,ylab=varname2,zlab="",
          main=paste0(varname1,":",varname2),theta=-30,phi=30,
          col="lightblue",lwd=.1)
  }
}
