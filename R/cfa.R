#' @title compute.cfa
#'
#' @description does the heavy lifting for computing cfa results
#' @inheritParams cfa
#'
#' @return CFA object
#' @keywords internal
#'
#' @export
compute.cfa <- function(tvals, yvals, data, yname, tname, xnames=NULL, drobj=NULL) {
    xmat <- data[,xnames]
    if (is.null(drobj)) {
        drobj <- TempleMetrics::distreg(yvals, data, yname,
                                        xnames=c(tname,xnames))
    }

    coef <- t(sapply(drobj$glmlist, coef))
    
    out <- list()
    
    for (i in 1:length(tvals)) {
        xmat1 <- cbind.data.frame(tvals[i], xmat)
        thisdist <- unlist(lapply(lapply(yvals, TempleMetrics::predict.DR, drobj=drobj, xdf=xmat1), mean))
        out[[i]] <- BMisc::makeDist(y_0, thisdist, rearrange=TRUE)##pred(tvals[i], drobj=drobj, yvals=yvals, xmat=xmat1)
    }

    return(CFA(tvals, out, coef=coef))

}

#' @title cfa
#'
#' @description compute counterfactuals using distribution regression
#'  with a continuous treatment
#'
#' @param tvals the values of the continuous treatment with which
#'  to compute counterfactuals
#' @param yvals the values to compute the counterfactual distribution fo
#' @param yname the name of the outcome (y) variable
#' @param tname the name of the treatment (t) variable
#' @param xnames the names of additional control variables to include
#' @param drobj optional distribution regression object that has been previously
#'  computed
#' @param se whether or not to compute standard errors using the bootstrap
#' @param iters how many bootstrap iterations to use
#' @param cl how many clusters to use for parallel computation of standard
#'  errors
#'
#' @return CFA object
#'
#' @export
cfa <- function(tvals, yvals, data, yname, tname, xnames=NULL, drobj=NULL,
                se=TRUE, iters=100, cl=1) {

    cfa.res <- compute.cfa(tvals, yvals, data, yname, tname, xnames, drobj)

    bootiterlist <- list()
    tvallist <- list()
    
    if (se) {
        cat("boostrapping standard errors...\n")
        n <- nrow(data)
        ##pb <- progress::progress_bar$new(total=iters)
        ##for (i in 1:iters) {
        ##    pb$tick()
        bstrap <- pbapply::pblapply(1:iters, function(z) {
            b <- sample(1:n, n, replace=TRUE)
            bdta <- data[b,]
            list(bootiter=compute.cfa(tvals, yvals, bdta, yname,
                                          tname, xnames, drobj),##$distcondt,
                 tvals=bdta[,tname])
        }, cl=cl)
        bootiterlist <- lapply(bstrap, function(x){ x$bootiter })
        tvallist <- lapply(bstrap, function(x) { x$tvals })
    }

    out <- CFA(tvals, cfa.res$distcondt, bootiterlist, tvallist, coef=cfa.res$coef)
}

#' @title title
#'
#' @description CFA objects
#'
#' @inheritParams cfa
#'
#' @param distcondt an ecdf object for a particular value of the treatment
#' @param bootiterlist a list of bootstrapped CFA objects that can be used
#'  for computing standard errors
#' @param tvallist the values of the treatment used in each bootstrap iteration
#'
#' @return CFA object
#' @export
CFA <- function(tvals, distcondt, bootiterlist=NULL, tvallist=NULL, coef=NULL) {
    out <- list(tvals=tvals, distcondt=distcondt, bootiterlist=bootiterlist, tvallist=tvallist, coef=coef)
    class(out) <- "CFA"
    out
}


#' @title getRes.CFA
#'
#' @description get a particular parameter of interest from a cfa object
#'
#' @param cfaobj a CFA object
#' @param fun a function to apply for every value of the treatment in the
#'  cfaobj
#' @param se whether or not to compute standard errors
#' @param setype whether to comput "pointwise" or "uniform" standard errors
#' @param ... can pass additional arguments to fun using this argument
#'
#' @return CFASE object
#' @export
getRes.CFA <- function(cfaobj, fun, se=T,  ...) {
    tvals <- cfaobj$tvals
    discondt <- cfaobj$distcondt
    bootiterlist <- lapply(cfaobj$bootiterlist, function(x) { x$distcondt })

    out <- t(simplify2array(lapply(discondt, fun, ...)))
    out <- if (nrow(out)==1) as.numeric(out) else out

    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist)) {
            bootout[[i]] <- t(simplify2array(lapply(bootiterlist[[i]], fun, ...)))
        }
        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        cb <- unlist(lapply(bootout, function(x) { max(abs((x-out)/ses)) } ))
        c <- quantile(cb, .95, type=1)
    }
    return(CFASE(tvals=tvals,est=out, se=ses, c=c))
}

#' @title getResDiff.CFA
#'
#' @description get a particular parameter of interest from a cfa object
#'
#' @param cfaobj a CFA object
#' @param fun a function to apply for every value of the treatment in the
#'  cfaobj
#' @param se whether or not to compute standard errors
#' @param setype whether to compute "pointwise" or "uniform" standard errors
#' @param ... can pass additional arguments to fun using this argument
#'
#' @return CFASE object
#' @export
getResDiff.CFA <- function(cfaobj1, cfaobj2, fun, se=T, ...) {
    tvals <- cfaobj1$tvals
    distcondt1 <- cfaobj1$distcondt
    bootiterlist1 <- lapply(cfaobj1$bootiterlist, function(x) { x$distcondt } )
    distcondt2 <- cfaobj2$distcondt
    bootiterlist2 <- lapply(cfaobj2$bootiterlist, function(x) { x$distcondt } )

    out <- t(simplify2array(lapply(distcondt1, fun, ...))) - t(simplify2array(lapply(distcondt2, fun, ...)))
    out <- if (nrow(out)==1) as.numeric(out) else out

    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist1)) {
            bootout[[i]] <- t(simplify2array(lapply(bootiterlist1[[i]], fun, ...))) - t(simplify2array(lapply(bootiterlist2[[i]], fun, ...)))
        }
        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        cb <- unlist(lapply(bootout, function(x) { max(abs((x-out)/ses)) } ))
        c <- quantile(cb, .95, type=1)
    }
    return(CFASE(tvals=tvals,est=out, se=ses, c=c))
}


#' @title getCoef.CFA
#'
#' @description get a particular parameter of interest from a cfa object
#'
#' @param cfaobj a CFA object
#' @param fun a function to apply for every value of the treatment in the
#'  cfaobj
#' @param se whether or not to compute standard errors
#' @param setype whether to comput "pointwise" or "uniform" standard errors
#' @param ... can pass additional arguments to fun using this argument
#'
#' @return CFASE object
#' @export
getCoef.CFA <- function(cfaobj, yvals, se=T,  ...) {
    
    coef <- cfaobj$coef
    bootiterlist <- cfaobj$bootiterlist

    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        ##for (i in 1:length(bootiterlist)) {
        ##    bootout[[i]] <- t(sapply(bootiterlist[[i]], function(x) { x$coef }))
        ##o}
        bootout <- lapply(bootiterlist, function(x) { x$coef })

        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        ##cb <- unlist(lapply(boot, function(x) { max(abs((x-coef)/ses)) } ))
        ##c <- quantile(cb, .95, type=1)
    }
    outlist <- lapply((1:ncol(coef)), function(i) {
        CFASE(tvals=yvals, est=coef[,i], se=ses[,i]) })
    return(outlist)
}


## #' 
## diff.cfa <- function(tvals, yvals, data, yname, tname,
##                      xnames1=NULL, drobj1=NULL,
##                      xnames2=NULL, drobj2=NULL,
##                      se=TRUE, iters=100, fun, ...) {

##     cfa1 <- compute.cfa(tvals, yvals, data, yname, tname, xnames1, drobj2)
##     cfa2 <- compute.cfa(tvals, yvals, data, yname, tname, xnames2, drobj2)

##     diffest <- getRes.CFA(cfa1, fun, se=F, ...)$est - getRes.CFA(cfa2, fun, se=F, ...)$est

##     ses <- NULL
##     if (se) {
##         bootiterlist <- list()
##         n <- nrow(data)
##         cat("bootstrapping standard errors...\n")
##         bootiterlist <- pbapply::pblapply(1:iters, function(z) {
##             b <- sample(1:n, n, T)
##             dtab <- data[b,]
##             cfa1b <- compute.cfa(tvals, yvals, dtab, yname, tname, xnames1, drobj2)
##             cfa2b <- compute.cfa(tvals, yvals, dtab, yname, tname, xnames2, drobj2)
##             diffestb <- as.matrix(getRes.CFA(cfa1b, fun, se=F, ...)$est - getRes.CFA(cfa2b, fun, se=F, ...)$est)
##         }, cl=8)

##         ses <- apply(simplify2array(bootiterlist), c(1,2), sd)
##         ses <- if (nrow(ses)==1 | ncol(ses)==1) as.numeric(ses) else ses
##     }
##     return(CFASE(tvals=tvals,est=diffest, se=ses))
## }


#' @title cfa2
#'
#' @description the same as cfa method except it computes two results at the
#'  same time which allows one to conduct inference on their difference
#'
#' @inheritParams cfa
#' @param xnames1 the first set of x variables
#' @param xnames2 the second set of x variables
#' @param drobj1 the first distribution regression object, optional
#' @param drobj2 the second distribution regression object, optional
#'
#' @return list of two CFA objects
#'
#' @export
cfa2 <- function(tvals, yvals, data, yname, tname,
                 xnames1=NULL, drobj1=NULL,
                 xnames2=NULL, drobj2=NULL,
                 se=TRUE, iters=100, cl=1) {

    cfa1 <- compute.cfa(tvals, yvals, data, yname, tname, xnames1, drobj2)
    cfa2 <- compute.cfa(tvals, yvals, data, yname, tname, xnames2, drobj2)


    bootiterlist1 <- list()
    bootiterlist2 <- list()
    tvallist <- list()
    
    if (se) {
        cat("boostrapping standard errors...\n")
        n <- nrow(data)
        ##pb <- progress::progress_bar$new(total=iters)
        ##for (i in 1:iters) {
        ##    pb$tick()
        bstrap <- pbapply::pblapply(1:iters, function(z) {
            b <- sample(1:n, n, replace=TRUE)
            bdta <- data[b,]
            list(bootiter1=compute.cfa(tvals, yvals, bdta, yname,
                                       tname, xnames1, drobj1),##$distcondt,
                 bootiter2=compute.cfa(tvals, yvals, bdta, yname,
                                       tname, xnames2, drobj2),
                 tvals=bdta[,tname])
        }, cl=cl)
        bootiterlist1 <- lapply(bstrap, function(x){ x$bootiter1 })
        bootiterlist2 <- lapply(bstrap, function(x){ x$bootiter2 })
        tvallist <- lapply(bstrap, function(x) { x$tvals })
    }

    out <- list(cfa1=CFA(tvals, cfa1$distcondt, bootiterlist1, tvallist, coef=cfa1$coef),
                cfa2=CFA(tvals, cfa2$distcondt, bootiterlist2, tvallist, coef=cfa2$coef))
}

#' @title test.CFA
#'
#' @description test if a counterfactual distribution is equal to its average
#'  for all values of the treatment
#'
#' @param cfaobj a CFA object
#' @param fun which function to use
#' @param allt all values of t in the dataset
#' @param se whether or not to compute standard errors
#' @param ... additional parameters for the function fun
#'
#' @return CFASE object
#'
#' @export
test.CFA <- function(cfaobj, fun, allt, se=T,  ...) {
    tvals <- cfaobj$tvals
    discondt <- cfaobj$distcondt
    bootiterlist <- lapply(cfaobj$bootiterlist, function(x) { x$distcondt } )
    tvallist <- cfaobj$tvallist

    out <- t(simplify2array(lapply(discondt, fun, ...)))
    fcondt <- approxfun(tvals, out)
    out <- out - mean(fcondt(allt), na.rm=TRUE) ## this drops some observations in the extreme tails
    out <- if (nrow(out)==1) as.numeric(out) else out

    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist)) {
            bout <- t(simplify2array(lapply(bootiterlist[[i]], fun, ...)))
            bfcondt <- approxfun(tvals, bout)
            bout <- bout - mean(bfcondt(tvallist[[i]]), na.rm=T)
            bootout[[i]] <- bout
        }
        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        cb <- unlist(lapply(bootout, function(x) { max(abs((x-out)/ses)) } ))
        c <- quantile(cb, .95, type=1)
    }
    return(CFASE(tvals=tvals,est=out, se=ses, c=c))
}


#' @title CFASE
#'
#' @description creates object of class CFASE
#'
#' @inheritParams cfa
#' @param est the estimate of the parameter at each value of t
#' @param c optional critical value for uniform inference
#'
#' @return CFASE object
#'
#' @export
CFASE <- function(tvals, est, se=NULL, c=NULL) {
    if (is.null(c)) {
        c <- 1.96
    }
    out <- list(tvals=tvals, est=est, se=se, c=c)
    class(out) <- "CFASE"
    out
}

#' @title ggplot2.CFA
#'
#' @description function for plotting results from counterfactual analysis
#'  using ggplot2
#'
#' @import ggplot2
#'
#' @param cfaseobj a CFASE object to plot
#' @param setype whether to plot pointwise, uniform, or both standard errors
#' @param ylim optional y limits on the plot
#' @param xlabel optional x axis labels
#' @param ylabel optional y axis labels
#' @param legend boolean for whether or not to plot a legend (tends to look
#'  better with this option set to FALSE)
#'
#' @return ggplot2 object
#'
#' @export
ggplot2.CFA <- function(cfaseobj, setype="pointwise", ylim=NULL,
                        xlabel=NULL, ylabel=NULL, legend=FALSE) {
    tvals <- cfaseobj$tvals
    est <- cfaseobj$est
    se <- cfaseobj$se
    c <- cfaseobj$c
    library(ggplot2)
    cmat1 <- cbind.data.frame(tvals, est)
    cmat1$which <- "est"
    cmat2 <- cbind.data.frame(tvals, se)
    cmat2$which <- "se"
    colnames(cmat2) <- colnames(cmat1)
    cmat <- rbind.data.frame(cmat1, cmat2)

    cmat <- tidyr::gather(cmat, key=k, value=v, -tvals, -which)

    cmat <- tidyr::spread(cmat, which, v)

    p <- ggplot(cmat, aes(x=tvals, y=est, group=k)) +
        geom_line(aes(color=k))
    if (setype == "both" | setype=="pointwise") {
            p <- p + geom_line(aes(y=est+1.96*se), lty=3) + geom_line(aes(y=est-1.96*se), lty=3)
    }
    if (setype=="both" | setype=="uniform") {
            p <-  p + geom_line(aes(y=est+c*se), lty=2) + geom_line(aes(y=est-c*se), lty=2)
    }        
    p <- p + theme_bw()

    if (!is.null(ylim)) {
        p <- p + scale_y_continuous(limits=ylim)
    }
    if (!is.null(xlab)) {
        p <- p + xlab(xlabel)
    }
    if (!is.null(ylab)) {
        p <- p + ylab(ylabel)
    }
    if (!legend) {
        p <- p + theme(legend.position="none")
    }
    p
}

