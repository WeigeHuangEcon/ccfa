#' Make counterfactual distribution functions

#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param ct A counterfactual treatment level
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles


#' @return The counterfactual distribution functions
#' @export

# to obtain the counterfactual distribution functions

makeCountDisFs.qr=function(Y,X,t,ct,s=NULL,e=NULL){
  countQuantsMat=getCountQuant.qr(Y,X,t,ct,s,e)
  df=rep(0,length(y))
  for(j in 1:length(df)){
    indicator=matrix(0,nrow = length(Y),ncol = s)
    for(i in 1:s){
      indicator[,i]=ifelse(countQuantsMat[,i]<sort(Y)[j],1,0)
    }
    df[j] =e+(1-2*e)/s*mean(rowSums(indicator))
  }
  df
}
