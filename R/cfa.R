compute.cfa <- function(tvals, yvals, data, yname, tname, xnames=NULL, drobj=NULL) {
    xmat <- data[,xnames]
    if (is.null(drobj)) {
        drobj <- TempleMetrics::distreg(yvals, data, yname,
                                        xnames=c(tname,xnames))
    }
    out <- list()
    
    for (i in 1:length(tvals)) {
        xmat1 <- cbind(tvals[i], xmat)
        thisdist <- unlist(lapply(lapply(yvals, predict.DR, drobj=drobj, xdf=xmat1), mean))
        out[[i]] <- BMisc::makeDist(y_0, thisdist, rearrange=TRUE)##pred(tvals[i], drobj=drobj, yvals=yvals, xmat=xmat1)
    }

    return(CFA(tvals, out))

}

cfa <- function(tvals, yvals, data, yname, tname, xnames=NULL, drobj=NULL,
                se=TRUE, iters=100) {

    cfa.res <- compute.cfa(tvals, yvals, data, yname, tname, xnames, drobj)

    bootiterlist <- list()
    
    if (se) {
        cat("boostrapping standard errors...\n")
        n <- nrow(data)
        pb <- progress::progress_bar$new(total=iters)
        for (i in 1:iters) {
            pb$tick()
            b <- sample(1:n, n, replace=TRUE)
            bdta <- data[b,]
            bootiterlist[[i]] <- compute.cfa(tvals, yvals, bdta, yname,
                                         tname, xnames, drobj)$distcondt
        }
    }

    out <- CFA(tvals, cfa.res$distcondt, bootiterlist)
}

CFA <- function(tvals, distcondt, bootiterlist=NULL) {
    out <- list(tvals=tvals, distcondt=distcondt, bootiterlist=bootiterlist)
    class(out) <- "CFA"
    out
}


getRes.CFA <- function(cfaobj, fun, se=T, ...) {
    tvals <- cfaobj$tvals
    discondt <- cfaobj$distcondt
    bootiterlist <- cfaobj$bootiterlist

    out <- t(simplify2array(lapply(discondt, fun, ...)))
    out <- if (nrow(out)==1) as.numeric(out) else out

    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist)) {
            bootout[[i]] <- t(simplify2array(lapply(bootiterlist[[i]], fun, ...)))
        }
        se <- apply(simplify2array(bootout), c(1,2), sd)
        se <- if (nrow(se)==1) as.numeric(se) else se
    }
    return(CFASE(tvals=tvals,est=out, se=se))
}

CFASE <- function(tvals, est, se=NULL) {
    out <- list(tvals=tvals, est=est, se=se)
    class(out) <- "CFASE"
    out
}

ggplot2.CFA <- function(cfaseobj) {
    tvals <- cfaseobj$tvals
    est <- cfaseobj$est
    se <- cfaseobj$se
    library(ggplot2)
    cmat1 <- cbind.data.frame(tvals, est)
    cmat1$which <- "est"
    cmat2 <- cbind.data.frame(tvals, se)
    cmat2$which <- "se"
    colnames(cmat2) <- colnames(cmat1)
    cmat <- rbind.data.frame(cmat1, cmat2)

    cmat <- tidyr::gather(cmat, key=k, value=v, -tvals, -which)

    cmat <- tidyr::spread(cmat, which, v)

    ggplot(cmat, aes(x=tvals, y=est, group=k)) +
        geom_line(aes(color=k)) +
        geom_line(aes(y=est+1.96*se), lty=2) +
        geom_line(aes(y=est-1.96*se), lty=2) +
        theme_bw()
    
}

