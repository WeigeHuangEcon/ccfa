
Var <- function(edf) {
    tau <- seq(.01,.99,.01)
    E <- mean(quantile(edf, tau, type=1))
    V <- mean( (quantile(edf, tau, type=1) - E)^2 )
    V
}

IQR <- function(edf, t1, t2) {
    quantile(edf, t1, type=1) - quantile(edf, t2, type=1)
}

pov <- function(edf, povline) {
    edf(povline)
}
