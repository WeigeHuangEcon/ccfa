#' @title E
#'
#' @description compute expectations from ecdf object
#'
#' @param edf an ecdf
#'
#' @return scalar expected value
#'
#'
#' @export
E <- function(edf) {
    tau <- seq(.01,.99,.01)
    E <- mean(quantile(edf, tau, type=1))
    E
}

#' @title Var
#'
#' @description compute variance from ecdf object
#'
#' @param edf an ecdf
#'
#' @return scalar variance
#'
#'
#' @export
Var <- function(edf) {
    tau <- seq(.01,.99,.01)
    E <- mean(quantile(edf, tau, type=1))
    V <- mean( (quantile(edf, tau, type=1) - E)^2 )
    V
}

#' @title IQR
#'
#' @description compute interquantile range from ecdf object
#'
#' @param edf an ecdf
#' @param t1 upper quantile
#' @param t2 lower quantile
#'
#' @return scalar interquantile range
#'
#'
#' @export
IQR <- function(edf, t1, t2) {
    quantile(edf, t1, type=1) - quantile(edf, t2, type=1)
}

#' @title pov
#'
#' @description compute fraction below the poverty line from ecdf
#'
#' @param edf an ecdf
#' @param povline the poverty line
#'
#' @return scalar fraction below poverty line
#'
#' @export
pov <- function(edf, povline) {
    edf(povline)
}

#' @title rich
#'
#' @description compute fraction of "rich"
#'
#' @param edf an ecdf
#' @param richline the cutoff for being rich
#'
#' @return scalar fraction that are "rich"
#'
#' @export
rich <- function(edf, richline) {
    1 - edf(richline)
}
