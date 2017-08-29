#' Get conditional quantiles for a treatment

#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles

#' @return The conditional quantiles
#' @export


# to obtain conditional quantiles estimators, we estimate S linear quantile regressions of Y on X

getCondQuants.qr=function(Y,X,t,s=NULL,e=NULL){
  if(is.null(s)){
    s=nrow(X)
  }
  if(is.null(e)){
    e=0.05
  }
  library(quantreg)
  taus=seq(e,1-e,length.out = s)
  rq(Y~t+X,tau = taus )
}
