stdcum<-function(klist){
#' Standardize a list of cumulants to have marginal zero mean and unit variance.
#'
#' @param klist a list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.
#' @return A list with the same components, with zero mean, unit marginal variance, and adjusted cumulants.  No adjustment is made for correlation.
#' @export
   for(ii in 1:2) for(jj in 1:2) for(kk in 1:2) for(ll in 1:2)
      klist$k4[ii,jj,kk,ll]<-klist$k4[ii,jj,kk,ll]/sqrt(klist$k2[ii,ii]*klist$k2[jj,jj]*klist$k2[kk,kk]*klist$k2[ll,ll])
   for(ii in 1:2) for(jj in 1:2) for(kk in 1:2)
      klist$k3[ii,jj,kk]<-klist$k3[ii,jj,kk]/sqrt(klist$k2[ii,ii]*klist$k2[jj,jj]*klist$k2[kk,kk])
   rho<-klist$k2[1,2]/sqrt(prod(diag(klist$k2)))
   klist$k2<-array(c(1,rho,rho,1),c(2,2))
   klist$k1<-rep(0,2)
   return(list(k1=klist$k1,k2=klist$k2,k3=klist$k3,k4=klist$k4))
}
