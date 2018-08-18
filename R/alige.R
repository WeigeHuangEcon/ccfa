#' @title localIGE.inner
#'
#' @description local intergenerational elasticities
#'
#' @param t conditional at a value T=t
#' @param Y outcome variable
#' @param T treatment variable
#' @param Xmat covariates
#' @param h bandwidth

#' @return lige
#' @export

# local ige at t
localIGE.inner <- function(t, Y, T, Xmat, h) {
  devtools::install_github("bcallaway11/TempleMetrics")
  dd <- llscm(t,Y,T,Xmat,h)
  xbar <- as.matrix(apply(Xmat, 2, mean))
  # m <- ncol(xmat)
  m <- ncol(Xmat)
  d <- dd[(m+1):length(dd)]
  t(xbar)%*%d
}

#' @title localIGE
#' @description localIGE
#' @param tvals a grid of values of treatment variable
#' @inheritParams localIGE.inner
#' @return
#' @export
# local ige at a grid of t's
localIGE <- function(tvals, Y, T, Xmat,h,cl=1) {
  pbsapply(tvals, localIGE.inner, Y=Y, T=T, Xmat=Xmat,cl=cl,h)
}


# wild bootstrap

#' @title sdF
#' @description using wild bootstrap to obtain standard deviation
#' @param B number of bootstrap iterations
#' @inheritParams localIGE
#' @return sd
#' @export
# a function used to compute standard divation of ALIGE
sdF=function(B,tvals, Y, T, Xmat,h){

  llscm_wb <- function(t, Y, T, Xmat, h) {
    dd <- llscm(t, Y, T, Xmat, h)
    X <- cbind(Xmat, (T - t)*Xmat)
    y <- as.matrix(Y)
    if (is.null(h)) {
      h <- 1.06*sd(T)*n^(-1/5) ## check that this is right
    }
    K <- diag(TempleMetrics::k(T-t, h=h, type="gaussian"))
    u=y-X%*%dd
    v=c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
    b_v=rbinom(nrow(dta), 1, prob=(sqrt(5)+1)/(2*sqrt(5)))
    b_v=ifelse(b_v==1,v[1],v[2])
    yhat=X%*%dd + u*b_v
    d <- solve(t(X)%*%K%*%X)%*%t(X)%*%K%*%yhat
    d
  }

  localIGE.inner_wb <- function(t, Y, T, Xmat, h) {
    dd <- llscm_wb(t,Y,T,Xmat,h)
    xbar <- as.matrix(apply(Xmat, 2, mean))
    m <- ncol(Xmat)
    d <- dd[(m+1):length(dd)]
    t(xbar)%*%d
  }

  localIGE_wb <- function(tvals, Y, T, Xmat,h) {
    seq=tvals
    for (i in 1:length(tvals)) {
      seq[i]=localIGE.inner_wb(tvals[i], Y, T, Xmat, h)
    }
    seq
  }
  library(foreach)
  library(doParallel)
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  registerDoParallel(no_cores)
  bmat=foreach(1:B,.combine = rbind)  %dopar%
    localIGE_wb(tvals, Y, T, Xmat,h)
  stopImplicitCluster()
  sd=apply(bmat, 2, sd)
  sd
}

#' @title Plot_ALIGE
#' @description plot ALIGE with 95\% confidence intervals
#' @param ALIGE
#' @param sd_ALIGE SD of ALIGE
#' @return Plot of ALIGE
#' @export
#'
Plot_ALIGE=function(ALIGE,sd_ALIGE){
  library(ggplot2)
  ALIGE_95=ALIGE+1.96*sd_ALIGE
  ALIGE_5=ALIGE-1.96*sd_ALIGE
  res=as.data.frame(cbind(ALIGE,sd_ALIGE))
  ggplot() +
    geom_line(data = res, aes(x = tvals, y = ALIGE), color = "black",lwd=1.5) +
    geom_point(aes(tvals,ALIGE),lwd=1.8,show.legend=F)+
    #geom_point(aes(tvals,ALIGE_95,color=2),lwd=1.4)+
    geom_line(data = res, aes(x = tvals, y = ALIGE_95), color = "black",lty=2,lwd=1.3) +
    geom_line(data = res, aes(x = tvals, y = ALIGE_5), color = "black",lty=2,lwd=1.3) +
    scale_y_continuous(limits=c(0,1)) +
    xlab("t") +
    ylab("ALIGE") +
    theme_bw() +
    theme(legend.title=element_blank())+
    theme(legend.background = element_blank())
}


