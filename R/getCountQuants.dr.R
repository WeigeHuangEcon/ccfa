#' Obtain the counterfactual quantiles
#' @param CondDiFs The conditional distributions functions
#' @param y0 The cutoff values
#' @param X The control variables
#' @param cts A list of counterfactual treatments
#' @param quants The quantiles to be obtained
#' @export
#' @return The counterfactual quantiles
getCountQuants.dr=function(CondDiFs,y0,X,cts,quants=NULL){
  if(is.null(quants)){
    quants=c(0.1,0.5,0.9)
  }
  re_CountDisFs.dr=re_makeCountDisFs.dr(CondDiFs,y0,X,cts)
  library(BMisc)
  quantile=matrix(0,nrow =length(quants), ncol = ncol(re_CountDisFs.dr))
  q=list()
  for(i in 1:length(quants)){
    for(j in 1:ncol(re_CountDisFs.dr)){
      q=cbind(q,makeDist(y0,re_CountDisFs.dr[,j]))
      quantile[,j]=quantile(q[[j]],probs = quants,type = 1)
    }
  }
  Average=colMeans(quantile)
  quantile=rbind(Average,quantile)
  rownames(quantile)=c('mean',round(quants,digits = 2))
  quantile
}
