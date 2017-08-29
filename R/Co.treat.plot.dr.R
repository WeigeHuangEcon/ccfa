
#' Plot coefficients on treatment variable
#' @param CondDiFs The conditional distributions functions
#' @param y0 The cutoff values
#' @export
#'
#'
Co.treat.plot.dr=function(CondDiFs,y0){
  coefficients_matrix=matrix(0,nrow =nrow(as.matrix(CondDiFs[[1]]$coefficients)),ncol =length(CondDiFs)  )
  for(i in 1: nrow(as.matrix(CondDiFs[[1]]$coefficients))){
    for(j in 1:length(CondDiFs)){
      coefficients_matrix[i,j]=( CondDiFs[[j]]$coefficients)[i]
    }
  }
  plot(y0,coefficients_matrix[2,],type = 'l',xlab = "",ylab = '')
}
