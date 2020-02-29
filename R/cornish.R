#' Calculate the Cornish-Fisher two-dimensional quantile approximation
#'
#' @param alpha a vector with two components giving a univariate upper tail probabiliity, and a bivariate upper tail probability, for wwhich the corresponding quanties are desired.
#' @param klist a list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.
#' @param nn the sample size.  If nn is not one, it is presumed that the cumulants of the sampling distribution under investigation behave in teh usual way; that is, third cumulants are k3/sqrt(nn), and fourth cumulants are k4/(nn.
#' @param xv As an alternative to alpha, one may enter ordinates associated with bivariate normal tail probabilities desired, and the tail probabilities are recalcuated from these.
#' @param prestd Logical variable indicating whether third and fourth order cumulants are already standardized to unit marginal variance.
#'
#' @return Approximate quantiles.
#' @export
cornish2<-function(alpha,klist,nn=1,xv=NULL,prestd=FALSE){
   cat("In cornish2 nn",nn,"\n")
   if(prestd){
      i4<-klist$k4
      i3<-klist$k3
      i2<-klist$k2
   }else{
      k2<-klist$k2
      i4<-array(NA,dim(klist$k4))
      i3<-array(NA,dim(klist$k3))
      i2<-array(NA,dim(k2))
      for(i in 1:2) for(j in 1:2) for(k in 1:2) for(l in 1:2) 
         i4[i,j,k,l]<-klist$k4[i,j,k,l]/sqrt(k2[i,i]*k2[j,j]*k2[k,k]*k2[l,l]) 
      for(i in 1:2) for(j in 1:2) for(k in 1:2) 
         i3[i,j,k]<-klist$k3[i,j,k]/sqrt( k2[i,i]* k2[j,j]* k2[k,k]) 
      for(i in 1:2) for(j in 1:2) i2[i,j]<-k2[i,j]/sqrt( k2[i,i]* k2[j,j])
   }
#  cat("In cornish2 alpha",alpha,"\n")
   if(is.null(xv)) xv<-fun.givex(alpha,i2[2,1])
   cat("In cornish2 xv",xv,"alpha",alpha,"\n")
   out<-.Fortran("bivcorn",
      x1=as.double(xv[1]), x2=as.double(xv[2]),
      x1p=as.double(0.0), x1pp=as.double(0.0),
      x2p=as.double(0.0), x2pp=as.double(0.0),
      rho=as.double(i2[2,1]),
      k111=as.double(i3[1,1,1]), k112=as.double(i3[1,1,2]),
      k122=as.double(i3[1,2,2]), k222=as.double(i3[2,2,2]),
      k1111=as.double(i4[1,1,1,1]), k1112=as.double(i4[1,1,1,2]),
      k1122=as.double(i4[1,1,2,2]), k1222=as.double(i4[1,2,2,2]),
      k2222=as.double(i4[2,2,2,2]),
      n=as.integer(nn),alpha=as.double(alpha),intflg=as.logical(FALSE),
      PACKAGE="bivcornish")
#   cat("raw   ",c(out$x1,out$x2),"\n")
#   cat("1st c ",c(out$x1p,out$x2p)/sqrt(nn),"\n")
#   cat("2nd c ",c(out$x1pp,out$x2pp)/(2*nn),"\n")
    return(list(raw=xv,cor=xv+c(out$x1p,out$x2p)/sqrt(nn)+c(out$x1pp,out$x2pp)/(2*nn)))
}
