#' Plot the counterfactual distribution functions before rearrangement
#' @param y0 The cutoff values
#' @param CountDisFs.dr The counterfactual distribution functions before rearrangement
#' @export

CountDisFs.dr.plot=function(y0,CountDisFs.dr){
  matplot(y0,CountDisFs.dr,type="l",ylab = '',xlab = '')
}
