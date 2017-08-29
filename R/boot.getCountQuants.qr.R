
#' Obtain the bootstrap version of quantiles

#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param cts A list of counterfactual treatments
#' @param quants The quantiles to be obtain
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles
#' @param B The number of bootstrap repetions
#' @return The bootstrap version of quantiles
#' @export

# bootstrap
boot.getCountQuants.qr <- function(B=NULL,Y,X,t,cts,quants=NULL,s=NULL,e=NULL) {
  if(is.null(quants)){
    quants=c(0.1,0.5,0.9)
  }
  if(is.null(B)){
    B=10
  }
  boot.results_list=list(list())
  for(b in 1:B){
    n <- nrow(X)
    bdata<-cbind(Y,t,X)
    bdata <- as.data.frame(bdata[sample(1:n, n, replace=T),])
    Y=bdata$Y
    X=as.matrix(bdata[,c(colnames(X))])
    t=bdata$t
    boot.result=getCountQuants.qr(Y,X,t,cts,quants,s,e)
    colnames(boot.result)<-round(cts,1)
    boot.results_list[[b]]= boot.result
  }
  boot.results_list
}
