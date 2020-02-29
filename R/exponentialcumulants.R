exponentialcumulants<-function(amat=rbind(c(1,1,0),c(0,1,1))){
#' Calculate cumulants for the linear compbination of independent exponential random variables, for exhibition of bivariate Cornish Fisher expansion.
#'
#' @param amat a vector with two rows and as many columns as there are independent exponentials
#'
#' @return A list with compoents k1, k2, k3, and k4, representing first through fourth comulant matrices respectively, having 1 through four dimensions, each dimension 2, respectively.
#' @export
#  browser()
   uc<-c(1,1,2,6)
   k0<-0
   k1<-array(0,rep(dim(amat)[1],1))
   k2<-array(0,rep(dim(amat)[1],2))
   k3<-array(0,rep(dim(amat)[1],3))
   k4<-array(0,rep(dim(amat)[1],4))
   for(k in seq(dim(amat)[2])) for(i in seq(dim(amat)[1])) {
      k1[i]<-k1[i]+uc[1]*amat[i,k]
      for(l in seq(dim(amat)[1])) {
         k2[i,l]<-k2[i,l]+uc[2]*amat[i,k]*amat[l,k]
         for (n in seq(dim(amat)[1])){
            k3[i,l,n]<-k3[i,l,n]+amat[i,k]*amat[l,k]*amat[n,k]*uc[3]
            for(j in seq(dim(amat)[1])){
               k4[i,l,n,j]<-k4[i,l,n,j]+amat[j,k]*amat[i,k]*amat[l,k]*amat[n,k]*uc[4]
            }
         }
      }
   }
   cumulantlist<-list( k0=k0, k1=k1, k2=k2, k3=k3, k4=k4,amat=amat)
   return(cumulantlist)
}
