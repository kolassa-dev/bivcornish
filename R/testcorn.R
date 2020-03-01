#' Calculate quantities involved in a Cornish-Fisher two-dimensional quantile approximation, for use in assessing relative accuracy of this approximation.
#'
#' @param alpha a vector with two components giving the first univariate upper target tail probability, and the second giving the bivariate upper tail probability.
#' @param xv a vector with two components giving the the ordinates associated with the quantiles; this may alternatively be a matrix with two rows and multiple columns, and the approximation will be evaluated at the pair in each column.  This argument takes precedence over alpha.
#' @param cumulants a function that generates a list of bivariate cumulant arrays, labeled k1 through k4, with 1 through 4 dimensions respectively, with each dimension 2.
#' @param amat a matrix transforming independent random variables into components of a dependent vector.
#' @param randf a function giving random observations, to be combined using amat below, to get a set of random vectors with two compoents.
#' @param nsamp The number of samples contributing to the Monte Carlo approximation to tail probabilities and quantiles.
#' @param nn the sample size.  If nn is not one, it is presumed that the cumulants of the sampling distribution under investigation behave in teh usual way; that is, third cumulants are k3/sqrt(nn), and fourth cumulants are k4/nn.
#'
#' @return A list with two components:
#'    An array with three dimensions, the first indicating quantities for the first dimension of the sample space, and its marginal distribution, and the second representing the second component and the joint distribution.  The second represents 1. Target upper probability, 2. The associated bivariate normal quantile, 3. The Edgeworth tail approximation for the normal quantile, 4. The Monte Carlo tail approximation for the normal quantile, 5. The Cornish-Fisher quantile, 6. The Monte Carlo tail approximation for the Cornish-Fisher quantile, 7. mean, 8. Standard deviations, 9. quantiles.
#'    The Monte Carlo sample drawn for comparison purposes.
#' @importFrom stats rexp
#' @export
testcorn<-function(alpha=c(.10,.05),xv=NULL,cumulants=exponentialcumulants,amat=rbind(c(1,1,0),c(0,1,1)),randf=rexp,nsamp=50000,nn=1){
   klist<-if(is.function(cumulants)) cumulants(amat=amat) else cumulants
   rho<-klist$k2[1,2]/sqrt(prod(diag(klist$k2)))
   slist<-stdcum(klist)
#  cat("rho",rho,"\n")
   if(length(dim(alpha))<2){
      alpha<-array(alpha,c(2,1))
   }
#  browser()
   if(is.null(xv)){
       cat("Recalculating xv")
       xv<-fun.givex(alpha,rho)
   }
   if(length(dim(xv))<2){
      xv<-array(xv,c(2,1))
   }
   ninter<-dim(klist$amat)[2]#Number of indepent observations to linearly combine.
# Draw random observations
   rawrandom<-array(randf(ninter*nn*nsamp), c(ninter,nn,nsamp))
# Transform using amat and take averages.
   transrandom<-apply(apply(rawrandom, c(2,3),"%*%",t(klist$amat)), c(1,3),"mean")
# Standardize to zero mean and unit marginal variance.
   xtest<-apply(apply(transrandom,2,"-",klist$k1),2,"/",sqrt(diag(klist$k2)))
   out<-array(NA,c(2,9,dim(xv)[2]))
   dimnames(out)<-list(NULL,c("alpha","Normal Quantile","Edgeworth Tail Approximation for Normal Quantile","Monte Carlo Tail for Normal Quantile","Cornish Fisher Quantile","Monte Carlo Tail for the CF Quantile","mean","sd","True Quantile"),NULL)
#  cat("Before big loop\n")
#  browser()
   for(jj in seq(dim(xv)[2])){
      out[,1,jj]<-apply(apply(apply(xtest,2,">",xv[,jj]),2,"cumprod"),1,"mean")
      cc<-cornish2(out[,1,jj],klist,nn=nn)
      out[,2,jj]<-cc$raw
      out[,3,jj]<-bivtail(out[,2,jj],slist,alreadystand=F)
      out[,4,jj]<-apply(apply(apply(xtest,2,">",out[,2,jj]),2,"cumprod"),1,"mean")
      out[,5,jj]<-cc$cor
      out[,6,jj]<-apply(apply(apply(xtest,2,">",out[,5,jj]),2,"cumprod"),1,"mean")
      out[,7,jj]<-klist$k1
      out[,8,jj]<-sqrt(diag(klist$k2))
      out[,9,jj]<-xv[,jj]
#     cat("target",out[,1,jj]," ")
#     cat("probs",c(out[1,1,jj],1-out[2,1,jj]/out[1,1,jj]),"\n")
#     out[1,9,jj]<-quantile(xtest[1,],probs=1-out[1,1,jj])
#     if(out[2,1,jj]<out[1,1,jj]){
#        out[2,9,jj]<-quantile(xtest[2,xtest[1,]>=out[1,9,jj]],
#           probs=1-out[2,1,jj]/out[1,1,jj])
#     }else{
#        cat("Bad quant at jj",jj,"out",out[,1,jj],"\n")
#     }
   }
   return(list(tails=out,xtest=xtest))
}
#testcorn(xv=makebigxv(10))#Make the figure for the poster
