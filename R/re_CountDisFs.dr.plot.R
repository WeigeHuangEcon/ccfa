#' Plot the counterfactual distribution functions after rearrangement
#' @param y0 The cutoff values
#' @param re_CountDisFs.dr The counterfactual distribution functions after rearrangement
#' @export
re_CountDisFs.dr.plot=function(y0,re_CountDisFs.dr){
  matplot(y0,re_CountDisFs.dr,type="l",ylab = '',xlab = '')
}
