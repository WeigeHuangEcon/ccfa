#' Obtain the standard errors matrix for quantiles by bootstrap using quantile regression approach

#' @param Y The outcome variable
#' @param X The control variables
#' @param B The number of bootstrap repetitions
#' @param t The treatment variable
#' @param cts A list of counterfactual treatments
#' @param quants The quantiles to be obtain
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles

#' @return The standard errors matrix for quantiles using quantile regression approach
#' @export
#'
getSDMat.qr=function(B=NULL,Y,X,t,cts,quants=NULL,s=NULL,e=NULL){
  if(is.null(quants)){
    quants=c(0.1,0.5,0.9)
  }
  if(is.null(B)){
    B=10
  }
  boot.results=boot.getCountQuants.qr(B=B,Y=Y,X=X,t=t,cts=cts,quants=quants,s=s,e=e)
  # obtain standard deviation
  sdmat <- apply(simplify2array(boot.results), c(1,2), sd)
  sdmat
}

