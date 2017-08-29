#' Obtain the standard errors matrix for quantiles by bootstrap
#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param y0 The cutoff values
#' @param cts A list of counterfactual treatments
#' @param B The number of bootstrap repetions
#' @param quants The quantiles to be obtained
#' @export
#' @return The standard errors matrix for quantiles
#'
getSDMat.dr=function(Y,X,t,y0,cts,B,quants=NULL){
  if(is.null(quants)){
    quants=c(0.1,0.5,0.9)
  }
  if(is.null(B)){
    B=10
  }
  CountQuants.dr_list=list(list())
  for(b in 1:B){
    n <- nrow(X)
    bdata<-cbind(Y,t,X)
    bdata <- as.data.frame(bdata[sample(1:n, n, replace=T),])
    Y=bdata$Y
    X=as.matrix(bdata[,c(colnames(X))])
    t=bdata$t

    CondDiFs=getCondDiFs.dr(Y=Y,X=X,t=t,y0=y0)
    # re_CountDisFs.dr=re_makeCountDisFs.dr(CondDiFs=CondDiFs,y0=y0,X=X,cts=cts)
    CountQuants.dr=getCountQuants.dr(CondDiFs=CondDiFs,y0=y0,X=X,cts=cts,quants=quants)
    CountQuants.dr_list[[b]]=CountQuants.dr
  }
  sdmat.dr <- apply(simplify2array(CountQuants.dr_list), c(1,2), sd)
  sdmat.dr
}
