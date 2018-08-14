
#' @title ll.scm
#'
#' @description local linear estimator of smoothing coefficient model
#'
#'
#'
#' @inheritParams cfa
#' @param condDist optional pre-estimated conditional distribution function
#'
#' @return CFA.OBJ
#' @keywords internal
#' @export

library(TempleMetrics)
library(pbapply)
setwd("/Users/weigehuang/Dropbox/LocalIGE/data")
dta <- read.csv("IGMdata.csv")

## local linear smooth varying coefficient model
ll.svc <- function(t, Y, T, Xmat, h) {
  X <- cbind(Xmat, (T - t)*Xmat)
  y <- as.matrix(Y)
  n <- length(T)
  if (is.null(h)) {
    h <- 1.06*sd(T)*n^(-1/5) ## check that this is right
  }
  K <- diag(TempleMetrics::k(T-t, h=h, type="gaussian"))
  ## kf=function (z, h , type = "gaussian")
  ## {
  ##   u <- z/h
  ##   if (type == "gaussian") {
  ##     dnorm(u) * (abs(u) <= 1)
  ##   }
  ##   else if (type == "epanechnikov") {
  ##     0.75 * (1 - u^2) * (abs(u) <= 1)
  ##   }
  ## }
  K<-diag(kf(z=T-t,h,type="gaussian"))
  dd <- solve(t(X)%*%K%*%X)%*%t(X)%*%K%*%y
  dd
}

localIGE.inner <- function(t, Y, T, Xmat, h) {
  dd <- ll.svc(t,Y,T,Xmat,h)
  xbar <- as.matrix(apply(Xmat, 2, mean))
  # m <- ncol(xmat)
  m <- ncol(Xmat)
  d <- dd[(m+1):length(dd)]
  t(xbar)%*%d
}

localIGE <- function(tvals, Y, T, Xmat,h,cl=1) {
  pbsapply(tvals, localIGE.inner, Y=Y, T=T, Xmat=Xmat,cl=cl,h)
}



# wild bootstrap

ll.svc_wb <- function(t, Y, Treat, Xmat, h) {
  X <- cbind(Xmat, (Treat - t)*Xmat)
  y <- as.matrix(Y)
  n <- length(Treat)
  if (is.null(h)) {
    h <- 1.06*sd(T)*n^(-1/5) ## check that this is right
  }
  #K <- diag(TempleMetrics::k(Treat-t, h=1, type="gaussian"))
  kf=function (z, h , type = "gaussian")
  {
    u <- z/h
    if (type == "gaussian") {
      dnorm(u) * (abs(u) <= 1)
    }
    else if (type == "epanechnikov") {
      0.75 * (1 - u^2) * (abs(u) <= 1)
    }
  }
  K<-diag(kf(z=Treat-t,h,type="gaussian"))
  dd <- solve(t(X)%*%K%*%X)%*%t(X)%*%K%*%y
  u=y-X%*%dd
  v=c(-(sqrt(5)-1)/2,(sqrt(5)+1)/2)
  #rbinom(v, 1, prob=c((sqrt(5)+1)/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5))))
  b_v=rbinom(nrow(dta), 1, prob=(sqrt(5)+1)/(2*sqrt(5)))
  b_v=ifelse(b_v==1,v[1],v[2])
  yhat=X%*%dd + u*b_v
  d <- solve(t(X)%*%K%*%X)%*%t(X)%*%K%*%yhat
  d
}

localIGE.inner_wb <- function(t, Y, Treat, Xmat, h) {
  dd <- ll.svc_wb(t,Y,Treat,Xmat,h)
  xbar <- as.matrix(apply(Xmat, 2, mean))
  # m <- ncol(xmat)
  m <- ncol(Xmat)
  d <- dd[(m+1):length(dd)]
  t(xbar)%*%d
}

localIGE_wb <- function(tvals, Y, Treat, Xmat,h) {
  pbsapply(tvals, localIGE.inner_wb, Y=Y, Treat=Treat, Xmat=Xmat,cl=8,h=h)
}

sdF=function(B,tvals, Y, Treat, Xmat,h){
  library(foreach)
  library(doParallel)
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  registerDoParallel(no_cores)
  bmat=foreach(1:B,.combine = list,
               .multicombine = T)  %dopar%
    localIGE_wb(tvals, Y, Treat, Xmat,h)
  stopImplicitCluster()
  sd=apply(simplify2array(bmat), 1, sd)
  sd
}

dta$yearborn <- dta$yearborn - 1970
tvals <- seq(log(20000),log(140000),length.out=20)
tvals <- seq(9.5,12,length.out=20)
xmat <- cbind(1, dta$headsex, dta$headnonwhite,dta$headlesshs,dta$headcol,dta$sex,dta$headveteran,dta$yearborn)
#xmat <- cbind(1, dta$headsex, dta$headnonwhite)
localIGE(tvals, dta$lcfincome, dta$lfincome, xmat,h)

x1 <- as.matrix(rep(1,nrow(dta)))
localIGE(tvals, dta$lcfincome, dta$lfincome, x1,h)


B=100
sd_LIGE=sdF(B,tvals, dta$lcfincome, dta$lfincome, x1,h)
LIGE=localIGE(tvals, dta$lcfincome, dta$lfincome, x1,h)

pdf("lige_wb.pdf")
plot(tvals,LIGE,type = "l",ylim = c(min(LIGE-1.96*sd_LIGE),max(LIGE+1.96*sd_LIGE)),xlab = "t")
lines(tvals,LIGE+1.96*sd_LIGE,lty=2)
lines(tvals,LIGE-1.96*sd_LIGE,lty=2)
dev.off()

sd_ALIGE=sdF(B,tvals, dta$lcfincome, dta$lfincome, xmat,h)
ALIGE=localIGE(tvals, dta$lcfincome, dta$lfincome, xmat,h)

pdf("alige_wb.pdf")
plot(tvals,ALIGE,type = "l",ylim = c(min(ALIGE-1.96*sd_ALIGE),max(ALIGE+1.96*sd_ALIGE)),xlab = "t")
lines(tvals,ALIGE+1.96*sd_ALIGE,lty=2)
lines(tvals,ALIGE-1.96*sd_ALIGE,lty=2)
dev.off()


plot(tvals,ALIGE,type = "l",ylim = c(-0.1,1),xlab = "t")
lines(tvals,ALIGE+1.96*sd_ALIGE,lty=2)
lines(tvals,ALIGE-1.96*sd_ALIGE,lty=2)
lines(tvals,LIGE,col=2)
lines(tvals,LIGE+1.96*sd_LIGE,lty=2,col=2)
lines(tvals,LIGE-1.96*sd_LIGE,lty=2,col=2)

# analytical standard deviations
sd_ALIGE_a=






  ## cross validation
  crossval.inner <- function(h,Y,T,xmat) {
    m <- ncol(xmat)
    mean(Y - pbsapply(1:length(Y), function(i) {
      coefs <- ll.svc(T[i], Y[-i], T[-i], xmat[-i,], h)
      coefs <- as.matrix(coefs[1:m]) ##drop the derivative ones
      xx <- as.matrix(xmat[i,])
      t(xx)%*%coefs
    }, cl=8))
  }

dta1 <- subset(dta, lfincome > quantile(dta$lfincome, .02) & lfincome < quantile(dta$lfincome, .98))
xmat <- cbind(1, dta1$headsex, dta1$headrace)
x1 <- as.matrix(rep(1,nrow(dta1)))
crossval.inner(1, dta1$lcfincome, dta1$lfincome, xmat) # 0.004400292


hh <- 1.06*sd(dta1$lfincome)*nrow(dta1)^(-1/5)
u <- 10*hh
l <- hh/10
o <- optimize(crossval.inner, lower=l, upper=u, Y=dta1$lcfincome, T=dta1$lfincome, xmat=xmat)
# objective 0.004400292 ; minimum 0.8547025


hval <- o$minimum ## need to choose higher upper bound next time # 0.8547025
hval

########################
## You could bootstrap standard errors for the difference between
## LIGE and ALIGE
## just call it twice
## not implemented
#########################
localIGE2 <- function() {
}


########################
### for getting various components of LIGE, not used
########################
componentsIGE <- function(tvals, Y, T, Xmat, h=NULL) {
  pbsapply(tvals, function(t) {
    dd <- ll.svc(t,Y,T,Xmat,h)
    xbar <- as.matrix(apply(Xmat, 2, mean))
    m <- ncol(Xmat)
    d <- dd[(m+1):length(dd)]
    xbar*d
  })
}

dta$yearborn <- dta$yearborn - 1970
tvals <- seq(log(20000),log(140000),length.out=20)
xmat <- cbind(1, dta$headsex, dta$headwhite, dta$headveteran,
              dta$sex, dta$yearborn, dta$headhs, dta$headcol)
ligeX <- localIGE(tvals, dta$lcfincome, dta$lfincome, xmat, h=.3)

## not sure if we want to do this or not
compLige <- componentsIGE(tvals, dta$lcfincome, dta$lfincome, xmat, h=.3)

x1 <- as.matrix(rep(1,nrow(dta)))
lige <- localIGE(tvals, dta$lcfincome, dta$lfincome, x1, h=.3)


## ggplot code
library(ggplot2)
cmat <- cbind.data.frame(tvals, ligeX, "ALIGE")
colnames(cmat) <- c("tvals", "lige", "adj")
cmat1 <- cbind.data.frame(tvals, lige, "LIGE")
colnames(cmat) <- colnames(cmat1) <- c("tvals", "lige", "adj")
cmat <- rbind.data.frame(cmat, cmat1)
cmat$adj <- as.factor(cmat$adj)
pdf("alige.pdf")
ggplot(cmat, aes(tvals, lige, group=adj)) +
  geom_line(aes(color=adj, linetype=adj), lwd=1.3) +
  geom_point(aes(color=adj), lwd=1.4) +
  scale_y_continuous(limits=c(0,1)) +
  xlab("Log of Parents' Income") +
  ylab("LIGE") +
  theme_bw() +
  theme(legend.title=element_blank())
dev.off()

diffLige <- lige - ligeX
