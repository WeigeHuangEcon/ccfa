#' Plot the counterfactual quantiles
#' @param cts A list of counterfactual treatments. 'cts' must be chosen from CountQuants.dr matrix
#' @param CountQuants.dr The counterfactual quantiles
#' @param quans The quantiles to plot, 'quans' must be chosen from CountQuants.dr matrix
#' @export

CountQuants.dr.plot=function(cts,CountQuants.dr,quans){
  matplot(cts,t(CountQuants.dr[quans,]),
          ylab = '',xlab = '',type = 'l',col = 1:length(quans))
}
