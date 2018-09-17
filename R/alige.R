
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
localIGE.inner <- function(formla,xformla=NULL, data,t, h=NULL) {
  formla=as.formula(formla)
  YT=model.frame(terms(formla,data=data),data=data)
  Y=YT[,1]
  T=YT[,2]
  n=length(Y)
  if (is.null(h)) {
    h <- 1.06*sd(T)*n^(-1/5) ## check that this is right
  }

  if(is.null(xformla)){
    dd <- TempleMetrics::llscm(formla,xformla,data,t,h)
    dd[2]
  } else {
    xformla=as.formula(xformla)
    dd <- TempleMetrics::llscm(formla,xformla,data,t,h)
    X=model.frame(terms(xformla,data=data),data=data)
    Xmat=as.matrix(X)
    xbar <- as.matrix(apply(Xmat, 2, mean))
    m <- ncol(Xmat)
    d <- dd[(m+2):length(dd)]
    t(xbar)%*%d
  }
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
#' formla=lcfincome~lfincome
#' xformla=~hs+col
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

localIGE <- function(formla, xformla=NULL, data,tvals,h=NULL,cl=1) {
  cl=cl
 out<- pbapply::pbsapply(tvals, localIGE.inner,formla=formla, xformla=xformla, data=data ,h=h,cl=cl)
 class(out)<-"localIGE"
 out
}

#' @title summary_localIGE
#'
#' @description prints a summary of a \code{localIGE} object
#'
#' @param object an \code{localIGE} object
#' @param ... extra arguments
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' formla=lcfincome~lfincome
#' xformla=~hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' cl=1
#' data=igm
#' object=localIGE(formla=formla, xformla=xformla, data=data,tvals=tvals,h=h,cl=cl)
#' summary_localIGE (object)
#' @export
#'

summary_localIGE <- function(object, ...) {
  citation()
  object
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
#' formla=lcfincome~lfincome
#' xformla=~hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' data=igm
#' B=7
#' sdF(B,formla=formla, xformla=xformla, data=data,tvals=tvals,h=h)
#' @export
# a function used to compute standard divation of ALIGE

sdF=function(B,formla, xformla, data,tvals,h,cl=1){

  sdF_t=function(formla, xformla, data,t,h){
    coef=TempleMetrics::llscm(formla,xformla,data,t,h)
    formla=as.formula(formla)
    YT=model.frame(terms(formla,data=data),data=data)
    Y=YT[,1]
    T=YT[,2]
    n=length(Y)
    if (is.null(h)) {
      h <- 1.06*sd(T)*n^(-1/5) ## check that this is right
    }
    if(is.null(xformla)){
      X=cbind(1,T-t)
      u=Y-X%*%coef
      v=c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
      b_v=rbinom(nrow(X), 1, prob=(sqrt(5)+1)/(2*sqrt(5)))
      b_v=ifelse(b_v==1,v[1],v[2])
      yhat=X%*%coef + u*b_v
      K <- diag(TempleMetrics::k(T-t, h=h, type="gaussian"))
      d <- solve(t(X)%*%K%*%X)%*%t(X)%*%K%*%yhat
      d[2]
    } else {
      xformla=as.formula(xformla)
      X=model.frame(terms(xformla,data=data),data=data)
      Xmat=cbind(1,X,X*(T-t))
      Xmat=as.matrix(Xmat)
      u=Y-Xmat%*%coef
      v=c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
      b_v=rbinom(nrow(X), 1, prob=(sqrt(5)+1)/(2*sqrt(5)))
      b_v=ifelse(b_v==1,v[1],v[2])
      yhat=Xmat%*%coef + u*b_v
      K <- diag(TempleMetrics::k(T-t, h=h, type="gaussian"))
      dd <- solve(t(Xmat)%*%K%*%Xmat)%*%t(Xmat)%*%K%*%yhat
      xbar <- as.matrix(apply(X, 2, mean))
      m <- ncol(X)
      d <- dd[(m+2):length(dd)]
      t(xbar)%*%d
    }
  }
  sdF_ts<- function(tvals) {
    seq=tvals
    for (i in 1:length(tvals)) {
      seq[i]=sdF_t(formla=formla, xformla=xformla, data=data,t=tvals[i],h=h)
    }
    seq
  }

  # library(foreach)
  # library(doParallel)

  ## chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  ## if (nzchar(chk) && chk == "TRUE") {
  ##   # use 2 cores in CRAN/Travis/AppVeyor
  ##   num_workers <- 2L
  ## } else {
  ##   # use all cores in devtools::test()
  ##   num_workers <- parallel::detectCores()
  ## }
  ## no_cores <- num_workers - 1

  ## ##no_cores <- parallel::detectCores() - 1
  ## no_cores <- cl
  ## cl<-parallel::makeCluster(no_cores)
  ## doParallel::registerDoParallel(cl)
  ## #doParallel::registerDoParallel(no_cores)
  ## bmat=foreach(1:B,.combine = rbind)  %dopar%
  ##   sdF_ts(tvals)
    ## stopImplicitCluster()

    
  bmat <- simplify2array(pblapply(1:B, function(b) { sdF_ts(tvals) }))

  sd=apply(bmat, 1, sd)
  class(sd)<-"sdF"
  sd
}



#' @title LLIGE
#' @description Computes local intergenerational elasticities and its Standard Deviations
#' @param B number of bootstrap iterations
#' @inheritParams localIGE
#' @return LLIGE AND its SD
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' formla=lcfincome~lfincome
#' xformla=~hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' data=igm
#' B=7
#' cl=1
#' LLIGE(B,formla=formla, xformla=xformla, data=data,tvals=tvals,h=h,cl=cl)
#' @export
#'

LLIGE=function(B,formla, xformla, data,tvals,h,cl){
 llige=localIGE(formla, xformla, data,tvals,h,cl=cl)
# names(llige)="llige"
 sd=sdF(B,formla, xformla, data,tvals,h)
# names(sd)="SD"
 out=list(llige=llige,SD=sd)
 class(out)<-"LLIGE"
 out
}

#' @title summary_LLIGE
#'
#' @description prints a summary of a \code{LLIGE} object
#' @param object an \code{LLIGE} object
#' @param ... extra arguments
#' @return summary of \code{LLIGE}
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' formla=lcfincome~lfincome
#' xformla=~hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' data=igm
#' B=7
#' cl=1
#' object=LLIGE(B,formla=formla, xformla=xformla, data=data,tvals=tvals,h=h,cl=cl)
#' summary_LLIGE(object)
#' @export

summary_LLIGE <- function(object, ...) {
  citation()
  object
}



#' @title Plot_ALIGE
#' @description plot ALIGE with 95\% confidence intervals
#' @param ALIGE lige from
#' @param sd_ALIGE SD of ALIGE
#' @param ylim ranges of values for y axis
#' @param xlab name of x axis
#' @param ylab name of y axis
#' @inheritParams localIGE
#' @return Plot of ALIGE
#' @examples
#' data(igm)
#' igm$hs=ifelse(igm$HEDUC=="HS",1,0)
#' igm$col=ifelse(igm$HEDUC=="COL",1,0)
#' formla=lcfincome~lfincome
#' xformla=~hs+col
#' tvals=seq(quantile(igm$lfincome,probs = 0.1),quantile(igm$lfincome,probs = 0.9),length.out = 10)
#' h=1.2
#' data=igm
#' cl=1
#' B=7
#' ALIGE=localIGE(formla=formla, xformla=xformla, data=data,tvals=tvals,h=h,cl=cl)
#' sd_ALIGE=sdF(B,formla=formla, xformla=xformla, data=data,tvals=tvals,h=h)
#' Plot_ALIGE(tvals,ALIGE,sd_ALIGE,xlab="t",ylab="ALIGE",ylim=c(0,1))
#' @export
#'
Plot_ALIGE=function(tvals,ALIGE,sd_ALIGE,xlab,ylab,ylim=c(0,1)){
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
    scale_y_continuous(limits=ylim) +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw() +
    theme(legend.title=element_blank())+
    theme(legend.background = element_blank())
  P
}


