#' Convert counterfactual distribution functions to quantiles

#' @param Y The outcome variable
#' @param X The control variables
#' @param t The treatment variable
#' @param ct A counterfactual treatment level
#' @param quants The quantiles to be obtained
#' @param s The number of linear quantile regressions
#' @param e The trimming level to avoid estimation of tail quantiles
#'
#' @return The counterfactual quantiles for a treatment
#' @details This function
#' @description Calculate the counterfactual quantiles for a treatment
#' @export

DisFs2Quants.qr=function(Y,X,t,ct,quants=NULL,s=NULL,e=NULL){
  if(is.null(quants)){
    quants=c(0.1,0.5,0.9)
  }
  fs=makeCountDisFs.qr(Y,X,t,ct,s,e)
  # library(BMisc)
  makeDist <- function(x, Fx, sorted=FALSE) {
    if (!sorted) {
      tmat <- cbind(x, Fx)
      tmat <- tmat[order(x),]
      x <- tmat[,1]
      Fx <- tmat[,2]
    }

    retF <- approxfun(x, Fx, method="constant",
                      yleft=0, yright=1, f=0, ties="ordered")
    class(retF) <- c("ecdf", "stepfun", class(retF))
    assign("nobs", length(x), envir = environment(retF))
    retF
  }

  quantiles=quantile(makeDist(sort(y),fs),probs=quants,type=1)
  mean=mean(quantiles)
  names(mean)='mean'
  quantiles=c(mean,quantiles)
  as.matrix(quantiles)
}
