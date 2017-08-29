#' Obtain the counterfacutal quantiles for one treatment
#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param ct A counterfactual treatment level
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles


#' @return The counterfacutal quantiles for one treatment
#' @export
#' @description Calculate the counterfactual quantiles


# to obtain the conditional counterfactual quantiles
# t is the treatment variable, ct is the counterfactual treatments. Both of them can be a set of varibles, such as t=c('lfincome','age')

getCountQuant.qr=function(Y,X,t,ct,s=NULL,e=NULL){
  if(is.null(s)){
    s=nrow(X)
  }
  if(is.null(e)){
    e=0.05
  }
  thisfit=getCondQuants.qr(Y,X,t,s,e)
  t=ct
  newdata=cbind(1,t,X)
  countQuantsMat=predict(thisfit,newdata=as.data.frame(newdata))
  countQuantsMat
}
