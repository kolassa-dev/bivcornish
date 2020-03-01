errorplot<-function(out,nn=1,mytitle="Exponential Example"){
#' Plot quantities representing errors in Cornish-Fisher approximation, relative to uncorrected normal theory approximation.
#'
#' @param out a matrix, generally provided by testcorn.  See testcorn documentation for details.
#' @param nn the sample size.  If nn is not one, it is presumed that the cumulants of the sampling distribution under investigation behave in teh usual way; that is, third cumulants are k3/sqrt(nn), and fourth cumulants are k4/(nn.
#' @param mytitle Title for plot.
#'
#' @return None.
#' @examples
#' #Generate Fig 1 of Chen, Kolassa, Seifu, and Zhong.  Next line takes 8 hours.  Uncommented example
#' #is not big enough to be useful, but runs quickly.
#' # errorplot(testcorn(xv=makebigxv(100),nsamp=500000,nn=1)$tail)
#' errorplot(testcorn(xv=makebigxv(10),nsamp=5000,nn=1)$tail)
#' @importFrom graphics abline contour legend lines par plot title
#' @importFrom stats quantile 
#' @export
   nps<-sqrt(dim(out)[3])
   firstdim<-out[1,,1:nps]
   par(mfrow=c(2,2))
   if(!is.null(mytitle)) par(oma=c(0,0,2,0)) else par(oma=c(0,0,0,0))
   rx<-range(firstdim[1,])
   ry<-range(abs(c(firstdim[1,]-firstdim[4,],firstdim[3,]-firstdim[4,])))
   plot(rx,ry,
      main="a. Absolute Error for Approximations",
      xlab=dimnames(firstdim)[[1]][1],ylab="Error", type="n")
   lines(firstdim[1,],abs(firstdim[1,]-firstdim[4,]),lty=1)
   lines(firstdim[1,],abs(firstdim[3,]-firstdim[4,]),lty=2)
   legend(rx[1],ry[2],lty=1:2,c("Normal","Edgeworth 4"))
   mm<-matrix(abs(out[2,1,]-out[2,4,])-abs(out[2,3,]-out[2,4,]),nrow=nps)
   xx<-out[1,2,seq(nps)]
   yy<-out[2,2,seq(nps)*nps]
   contour(xx,yy,mm,xlab="First Ordinate",ylab="Second Ordinate",
      main="b. Difference in Absolute Error", sub="Normal Minus Edgeworth")
   contour(xx,yy,matrix(out[2,4,],nrow=nps),
      xlab="First Ordinate",ylab="Second Ordinate",
      levels=c(.5,.25,.1,.05,.01,.005,.001,.0005,.0001),
      main="c. True Bivariate Probability")
   v<-unique(out[2,9,])
   uln<-ulc<-umn<-umc<-rep(NA,length(v))
   for(j in seq(length(v))){
      uln[j]<-quantile(out[2,2,out[2,9,]==v[j]],.01)
      ulc[j]<-quantile(out[2,5,out[2,9,]==v[j]],.01)
      umn[j]<-quantile(out[2,2,out[2,9,]==v[j]],.99)
      umc[j]<-quantile(out[2,5,out[2,9,]==v[j]],.99)
   }
   plot(range(v),range(c(uln,ulc,umn,umc)),type="n",
      xlab="True Ordinate",ylab="Approximate Ordinate",
      main="d. Second Ordinate Approximation\n Range")
   lines(v,uln,lty=2,col=2);lines(v,umn,lty=2,col=2)
   lines(v,ulc,lty=3,col=3);lines(v,umc,lty=3,col=3)
   abline(a=0,b=1)
   legend(min(v),max(c(umn,umc)),lty=c(2,0,3),col=c(2,0,3),
      legend=c("Multivariate","Normal","New"))
   if(!is.null(mytitle)) title(paste(mytitle,"sample size",nn),outer=TRUE)
}
