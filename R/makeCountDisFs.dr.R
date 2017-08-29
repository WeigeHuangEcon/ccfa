#' Obtain the counterfactual distribution functions before rearrangement
#' @param CondDiFs The conditional distributions functions
#' @param y0 The cutoff values
#' @param X The control variables
#' @param cts A list of counterfactual treatments
#' @export
#' @return The counterfactual distribution functions before rearrangement
#'
makeCountDisFs.dr=function(CondDiFs,y0,X,cts){
  CountDisFs.dr=matrix(0,nrow = length(y0),ncol = length(cts))
  for (i in 1: length(y0)){
    for(j in 1:length(cts)){
      predictinput=data.frame(cbind(cts[j],X))
      names(predictinput)<-c(names(CondDiFs[[i]]$model)[2],colnames(X))
      CountDisFs.dr[i,j]=mean(predict(CondDiFs[[i]],predictinput,type = 'response'))
    }
  }
  CountDisFs.dr
}
