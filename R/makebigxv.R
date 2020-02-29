#' Calculate a grid of points on which to test a Cornish-Fisher expansion.
#'
#' @param nps number of points (less one) on for each dimension of the resulting grid.
#'
#' @return Grid of ordinates.
#' @export
makebigxv<-function(nps){
   xxx<-0:nps
   return(2*rbind(rep(xxx,nps+1),rep(xxx,rep(nps+1,nps+1)))/nps)
}
