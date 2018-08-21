
#' @title localIGE.inner
#'
#' @description local intergenerational elasticities
#' @param formla a formula y ~ treatment
#' @param xformla one sided formula for x variables to include, e.g. ~x1 + x2
#' @param data the data.frame where y, t, and x are
#' @param t conditional at a value T=t (i.e. should be scalar)
#' @param h bandwidth
#' @return lige
#'
#' @keywords internal
#' @export
localIGE.inner <- function(formla,xformla, data,t, h) {
  dd <- TempleMetrics::llscm(formla,xformla,data,t,h)
  X=model.frame(terms(xformla,data=data),data=data)
  Xmat=as.matrix(X)
  xbar <- as.matrix(apply(Xmat, 2, mean))
  m <- ncol(Xmat)
  d <- dd[(m+1):length(dd)]
  t(xbar)%*%d
}

#' @title localIGE
#' @description Computes local intergenerational elasticities
#' @param tvals a grid of values of treatment variable
#' @param cl the number of clusters to use, default is 1
#' @inheritParams localIGE.inner
#' @return lige
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' igm$const=1
#' formla=lcfincome~lfincome
#' xformla=~const+hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' cl=1
#' data=igm
#' localIGE(formla=formla, xformla=xformla, data=data,tvals=tvals,h=h,cl=cl)
#' @export


# local ige at a grid of t's
#localIGE <- function(tvals, Y, T, Xmat,h,cl=1) {
#  pbapply::pbsapply(tvals, localIGE.inner, Y=Y, T=T, Xmat=Xmat,cl=cl,h)
#}

localIGE <- function(formla, xformla, data,tvals,h,cl=1) {
  formla=as.formula(formla)
  xformla=as.formula(xformla)
  YT=model.frame(terms(formla,data=data),data=data)
  X=model.frame(terms(xformla,data=data),data=data)
  Y=YT[,1]
  T=YT[,2]
  Xmat=as.matrix(X)
  cl=cl
  pbapply::pbsapply(tvals, localIGE.inner,formla=formla, xformla=xformla, data=data ,h=h,cl=cl)
}

# wild bootstrap

#' @title sdF
#' @description using wild bootstrap to obtain standard deviation
#' @param B number of bootstrap iterations
#' @inheritParams localIGE
#' @return sd
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' igm$const=1
#' formla=lcfincome~lfincome
#' xformla=~const+hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' data=igm
#' B=7
#' sdF(B,formla=formla, xformla=xformla, data=data,tvals=tvals,h=h)
#' @export
# a function used to compute standard divation of ALIGE
sdF=function(B,formla, xformla, data,tvals,h){
  formla=as.formula(formla)
  xformla=as.formula(xformla)
  YT=model.frame(terms(formla,data=data),data=data)
  X=model.frame(terms(xformla,data=data),data=data)
  Y=YT[,1]
  T=YT[,2]
  Xmat=as.matrix(X)

  llscm_wb <- function(t, Y, T, Xmat, h) {
    dd <- TempleMetrics::llscm.inner(t, Y, T, Xmat, h)
    X <- cbind(Xmat, (T - t)*Xmat)
    y <- as.matrix(Y)
    if (is.null(h)) {
      h <- 1.06*sd(T)*n^(-1/5) ## check that this is right
    }
    K <- diag(TempleMetrics::k(T-t, h=h, type="gaussian"))
    u=y-X%*%dd
    v=c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
    b_v=rbinom(nrow(X), 1, prob=(sqrt(5)+1)/(2*sqrt(5)))
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
  #library(foreach)
  #library(doParallel)
  no_cores <- parallel::detectCores() - 1
  cl<-parallel::makeCluster(no_cores)
  doParallel::registerDoParallel(cl)
  #doParallel::registerDoParallel(no_cores)
  bmat=foreach(1:B,.combine = rbind)  %dopar%
    localIGE_wb(tvals, Y, T, Xmat,h)
  stopImplicitCluster()

  sd=apply(bmat, 2, sd)
  sd
}

#' @title Plot_ALIGE
#' @description plot ALIGE with 95\% confidence intervals
#' @param ALIGE lige from
#' @param sd_ALIGE SD of ALIGE
#' @param xlab name of x axis
#' @param ylab name of y axis
#' @return Plot of ALIGE
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' igm$const=1
#' formla=lcfincome~lfincome
#' xformla=~const+hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' data=igm
#' cl=1
#' B=7
#' ALIGE=localIGE(formla=formla, xformla=xformla, data=data,tvals=tvals,h=h,cl=cl)
#' sd_ALIGE=sdF(B,formla=formla, xformla=xformla, data=data,tvals=tvals,h=h)
#' Plot_ALIGE(ALIGE,sd_ALIGE,xlab="t",ylab="ALIGE")
#' @export
#'
Plot_ALIGE=function(ALIGE,sd_ALIGE,xlab,ylab){
  #library(ggplot2)
  ALIGE_95=ALIGE+1.96*sd_ALIGE
  ALIGE_5=ALIGE-1.96*sd_ALIGE
  res=as.data.frame(cbind(ALIGE,sd_ALIGE))
  P=ggplot() +
    geom_line(data = res, aes(x = tvals, y = ALIGE), color = "black",lwd=1.5) +
    geom_point(aes(tvals,ALIGE),lwd=1.8,show.legend=F)+
    #geom_point(aes(tvals,ALIGE_95,color=2),lwd=1.4)+
    geom_line(data = res, aes(x = tvals, y = ALIGE_95), color = "black",lty=2,lwd=1.3) +
    geom_line(data = res, aes(x = tvals, y = ALIGE_5), color = "black",lty=2,lwd=1.3) +
    scale_y_continuous(limits=c(0,1)) +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw() +
    theme(legend.title=element_blank())+
    theme(legend.background = element_blank())
  P
}


