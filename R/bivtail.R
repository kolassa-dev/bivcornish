#' Calculate the Edgworth two-dimensional tail probability approximation
#'
#' @param xv a vector with two components giving the point at which to evaluate the tail approximaton
#' @param klist a list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.
#' @param nn the sample size.  If nn is not one, it is presumed that the cumulants of the sampling distribution under investigation behave in teh usual way; that is, third cumulants are k3/sqrt(nn), and fourth cumulants are k4/(nn.
#' @param alreadystand Logical flag, false if ordinates need to be standardized to zero mean and unit variance, and true otherwise.
#'
#' @return A two-component vector with the marginal univariate Edgeworth approximation to the tail for the first dimension, and the bivariate tail approximation.
#' @importFrom mvtnorm pmvnorm
#' @export
bivtail<-function(xv,klist,nn=1,alreadystand=F){
   if(!alreadystand) xv<-(xv-klist$k1)/sqrt(diag(klist$k2))
   klist<-stdcum(klist)
   ttt<-c(0,pmvnorm(lower=as.vector(xv),corr=klist$k2))
   dd<-.Fortran("bivtail",
      x1=as.double(xv[1]), x2=as.double(xv[2]),
      tail=as.double(ttt),rho=as.double(klist$k2[1,2]), 
      k111=as.double(klist$k3[1,1,1]), k112=as.double(klist$k3[1,1,2]),
      k122=as.double(klist$k3[1,2,2]), k222=as.double(klist$k3[2,2,2]),
      k1111=as.double(klist$k4[1,1,1,1]), k1112=as.double(klist$k4[1,1,1,2]),
      k1122=as.double(klist$k4[1,1,2,2]), k1222=as.double(klist$k4[1,2,2,2]), 
      k2222=as.double(klist$k4[2,2,2,2]),nn=as.integer(nn),PACKAGE="bivcornish")
   return(dd$tail)
}
