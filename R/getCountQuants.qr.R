#' Get conditional quantiles for a list of treatments

#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param cts A list of counterfactual treatments
#' @param quants The quantiles to be obtain
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles



#' @return The conditional quantiles for a list of treatments
#' @export


# obtain unconditional counterfactual quantiles

getCountQuants.qr <- function(Y,X,t,cts,quants=NULL,s=NULL,e=NULL) {
  if(is.null(quants)){
    quants=c(0.1,0.5,0.9)
  }
  n_treatment=length(cts)
  cQMat=DisFs2Quants.qr(Y,X,t,ct=cts[1],quants,s,e)
  for(n in 2:n_treatment){
    m=DisFs2Quants.qr(Y,X,t,ct=cts[n],quants,s,e)
    cQMat=cbind(cQMat,m)
  }
  colnames(cQMat)<-round(cts,1)
  cQMat
}
