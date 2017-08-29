#' Obtain the conditional distribution functions
#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param y0 A vector of cutoff values
#' @export
#'
#' @return The conditional distribution functions
#'
# obtain conditional distribution functions
getCondDiFs.dr=function(Y,X,t,y0){
  # create the indicator dependent variables '1(Y <= y_0)' for each 'y0'
  indicatormatrix=matrix(0,nrow = nrow(X),ncol =length(y0) )
  for(i in 1:nrow(X)){
    for(j in 1:length(y0))
      if (Y[i] <= y0[j]) {
        indicatormatrix[i,j]=1
      }
    else indicatormatrix[i,j]=0
  }

  CondDiFs=list(list())
  for(i in 1:ncol(indicatormatrix)){
    CondDiFs[[i]] <- glm(indicatormatrix[,i] ~t+x,family=binomial(link='logit'))
  }
  CondDiFs
}
