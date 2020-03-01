#' Give the bivariate normal upper quantile
#'
#' @param alpha a vector with two components giving upper quantiles
#' @param rho the correlation for the bivariate normal with two normal components.
#' @param expect the two-component vector of expected values of the bivariate distribution, if not zero.
#' @param sd the two-component vector of marginal standard deviations of the bivariate distribution, if not 1.
#'
#' @return A two-component vector whose first component is the upper standard normal quantile associated with alpha[1], and whose second component is that one to make the bivariate normal upper tail area alpha[2].
#' @importFrom stats qnorm
#' @export
fun.givex<-function(alpha,rho,expect=NULL,sd=NULL){
#  cat("alpha",alpha,"rho",rho,"\n")
   x1<-qnorm(1-alpha[1])
   out<-.Fortran("isecnorm",x1=as.double(x1),x2=as.double(0),targ=as.double(alpha[2]),rho=as.double(rho),efg=as.integer(0),PACKAGE="bivcornish")
   if(!is.null(expect)){
      out$newx<-expect+sd*c(out$x1,out$x2)
   }
   return(c(out$x1,out$x2))
}
