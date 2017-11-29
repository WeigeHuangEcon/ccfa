#' @title compute.cfa2
#'
#' @description an update of compute.cfa that allows one to use other
#'  other estimators such as quantile regression or nonparametric methods
#'  (not yet implemented) as the first step estimator of the conditional
#'  distribution
#'
#' @inheritParams cfa
#' @param condDist optional pre-estimated conditional distribution function
#'
#' @return CFA.OBJ
#' @keywords internal
#' @export
compute.cfa2 <- function(tvals, yvals, data, yname, tname, xnames=NULL,  method="dr", link="logit", tau=seq(.1,.99,.01), condDistobj=NULL) {
    obj <- condDistobj
    xmat <- data[,xnames]
    if (is.null(obj)) {
        formla <- as.formula(paste0(yname,"~",tname))
        formla <- addCovToFormla(xnames, formla)
        if (method == "dr") {
            obj <- distreg(formla, data, yvals, link)
        } else if (method == "qr") {
            obj <- rq(formla, tau=tau, data)
        } else if (method == "ll") {
            formla <- as.formula(paste0(yname,"~",tname))
            xformla <- as.formula("y ~ xxxx")
            xformla <- addCovToFormla(xnames, xformla)
            xformla <- dropCovFromFormla("xxxx", xformla)
            formula.tools::lhs(xformla) <- NULL
            print("computing local linear distribution regression")
            obj <- lldistreg(formla, xformla, data, yvals, tvals)
        } else {
            stop("method not yet implemented")
        }
    }

    coef <- NULL
    ##coef <- t(sapply(drobj$glmlist, coef))
    
    out <- list()
    
    ##for (i in 1:length(tvals)) {

    out <- lapply(1:length(tvals), function(i) { 
        xmat1 <- data[,c(tname,xnames)]
        xmat1 <- cbind.data.frame(1, xmat1)
        xmat1[,tname] <- tvals[i]
        thisdist <- Fycondx(obj, yvals, xmat1)
        combineDfs(yvals, thisdist)
    } )

    return(CFA.OBJ(tvals, out, coef=coef))

}


#' @title cfa.inner
#'
#' @description calls function to compute counterfactuals
#'
#' @param yname the name of the outcome (y) variable
#' @param tname the name of the treatment (t) variable
#' @param xnames the names of additional control variables to include
#' @inheritParams cfa
#'
#' @return CFA object
#'
#' @keywords internal
#' 
#' @export
cfa.inner <- function(tvals, yvals, data, yname, tname, xnames=NULL,
                      method="dr", link="logit", tau=seq(.01,.99,.01),
                      condDistobj=NULL,
                      se=TRUE, iters=100, cl=1) {

    cfa.res <- compute.cfa2(tvals, yvals, data, yname, tname, xnames,
                            method, link, tau, condDistobj)

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
            list(bootiter=compute.cfa2(tvals, yvals, bdta, yname,
                                      tname, xnames, method,
                                      link, tau),##$distcondt,
                 tvals=bdta[,tname])
        }, cl=cl)
        bootiterlist <- lapply(bstrap, function(x){ x$bootiter })
        tvallist <- lapply(bstrap, function(x) { x$tvals })
    }

    out <- CFA.OBJ(tvals, cfa.res$distcondt, bootiterlist, tvallist, coef=cfa.res$coef)
}

#' @title cfa
#'
#' @description compute counterfactuals using distribution regression
#'  with a continuous treatment
#'
#' @param formla a formula y ~ treatment
#' @param xformla one sided formula for x variables to include, e.g. ~x1 + x2
#' @param tvals the values of the "treatment" to compute parameters of
#'  interest for
#' @param yvals the values to compute the counterfactual distribution for
#' @param data the data.frame where y, t, and x are
#' @param method either "dr" or "qr" for distribution regression or quantile regression
#' @param link if using distribution regression, any link function that works with the binomial family (e.g. logit (the default), probit, cloglog)
#' @param tau if using quantile regression, which values of tau to estimate
#'  the conditional quantiles
#' @param condDistobj optional conditional distribution  object that has
#'  been previously  computed
#' @param se whether or not to compute standard errors using the bootstrap
#' @param iters how many bootstrap iterations to use
#' @param cl how many clusters to use for parallel computation of standard
#'  errors
#' 
#' @return CFA object
#'
#' @examples
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' ## This line doesn't adjust for any covariates
#' cfa(lcfincome ~ lfincome, tvals=tvals, yvals=yvals, data=igm,
#'  se=FALSE)
#'
#' ## This line adjusts for differences in education
#' cfa(lcfincome ~ lfincome, ~HEDUC, tvals=tvals, yvals=yvals, data=igm,
#'  se=FALSE)
#' 
#' @export
cfa <- function(formla, xformla=NULL, tvals, yvals, data,
                method="dr", link="logit", tau=seq(.01,.99,.01),
                condDistobj=NULL,
                se=TRUE, iters=100, cl=1) {

    formla <- as.formula(formla)
    dta <- model.frame(terms(formla,data=data),data=data) #or model.matrix
    yname <- colnames(dta)[1]
    tname <- colnames(dta)[2]
    xnames <- NULL
    ##set up the x variables
    if (!(is.null(xformla))) {
        xformla <- as.formula(xformla)
        xformla <- addCovToFormla("0", xformla)
        xdta <- model.matrix(xformla, data=data)
        dta <- cbind.data.frame(dta, xdta)
        xnames <- colnames(xdta)[-1]
    }

    cfa.inner(tvals, yvals, dta, yname, tname, xnames, method, link, tau,
              condDistobj, se, iters, cl)
}
    

#' @title CFA.OBJ
#'
#' @description CFA objects
#'
#' @inheritParams cfa
#'
#' @param distcondt an ecdf object for a particular value of the treatment
#' @param bootiterlist a list of bootstrapped CFA objects that can be used
#'  for computing standard errors
#' @param tvallist the values of the treatment used in each bootstrap iteration
#' @param coef the coefficients from a distribution regression
#'
#' @return CFA object
#' @export
CFA.OBJ <- function(tvals, distcondt, bootiterlist=NULL, tvallist=NULL, coef=NULL) {
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
#'  cfaobj.  The ccfa package provides several built-in functions:
#'  E (for expected value as a function of the treatment variable),
#'  Var (for the variance as a function of the treatment variable),
#'  IQR (the interquantile range as a function of the treatment variable),
#'  pov (the fraction of observations with outcomes below some threshold,
#'   as a function of the treatment variable),
#'  rich (the fraction of observations with outcomes above some threshold,
#'   as a function of the treatment variable),
#' but other user-defined functions can be written.  The requirement is that
#'  they need to take in an ecdf object and output a scalar result.
#' @param se whether or not to compute standard errors
#' @param ... can pass additional arguments to fun using this argument
#'
#' @return CFASE object
#' @examples
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95),
#'   length.out=50)
#' 
#' ## obtain counterfactual results
#' cfaresults <- cfa(lcfincome ~ lfincome, tvals=tvals, yvals=yvals, data=igm,
#'  se=FALSE)
#'
#' ## get the average outcome (lfincome) as a function of the treatment
#' ## variable (lfincome)
#' getRes.CFA(cfaresults, E, se=FALSE)
#'
#' ## get the variance of the outcomes as a function of the treatment
#' ## variable
#' getRes.CFA(cfaresults, Var, se=FALSE)
#'
#' ## get the inter-quantile range of outcomes as a function of the
#' ## treatment variable
#' getRes.CFA(cfaresults, IQR, se=FALSE, t1=0.9, t2=0.1)
#' 
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
#' @description Get the difference between two CFA objects
#'
#' @param cfaobj1 the first CFA object
#' @param cfaobj2 the second CFA object
#' @param fun a function to apply for every value of the treatment in the
#'  cfaobj
#' @param se whether or not to compute standard errors
#' @param ... can pass additional arguments to fun using this argument
#'
#' @return CFASE object
#' @examples
#' \dontrun{
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' 
#' ## obtain counterfactual results
#' out <- cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr",
#' xformla2=~HEDUC, method2="qr", iters=10, tau1=seq(.05,.95,.05),
#' tau2=seq(.05,.95,.05))
#'
#' ## get the difference between the average that adjusts for covariates and
#' ## the one that does not
#' getResDiff.CFA(out$cfa1, out$cfa2, E)
#' }
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
#' @param yvals the y values that the cfa object is computed for
#' @param se whether or not to compute standard errors
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




#' @title cfa2
#'
#' @description the same as cfa method except it computes two results at the
#'  same time which allows one to conduct inference on their difference
#'
#' @inheritParams cfa
#' @param xformla1 an optional formula for the first set of x variables
#' @param method1 the first method for estimating the conditional distribution
#'  it can be "dr" for distribution regression or "qr" for quantile regression
#' @param link1 if using distribution regression, set the link variable.  It
#'  can be any link function accepted by glm, e.g. logit, probit, cloglog
#' @param tau1 if using quantile regression, the values of tau to use, the
#'  default is seq(.01,.99,.01)
#' @param condDistobj1 if have already calculated a conditional distribution
#'  object outside of the model, can set it here
#' @param xformla2 an optional formula for the second set of x variables
#' @param method2 the second method for estimating the conditional distribution
#'  it can be "dr" for distribution regression or "qr" for quantile regression
#' @param link2 if using distribution regression, set the link variable.  It
#'  can be any link function accepted by glm, e.g. logit, probit, cloglog
#' @param tau2 if using quantile regression, the values of tau to use, the
#'  default is seq(.01,.99,.01)
#' @param condDistobj2 if have already calculated a conditional distribution
#'  object outside of the model, can set it here
#'
#' @return list of two CFA objects
#'
#' @examples
#' #' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' 
#' ## obtain counterfactual results using quantile regression with
#' ## no covariates and adjusting for education
#' cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr", xformla2=~HEDUC,
#' method2="qr", se=FALSE, tau1=seq(.1,.9,.1), tau2=seq(.1,.9,.1))
#'
#' @export
cfa2 <- function(formla, tvals, yvals, data, 
                 xformla1=NULL, method1="dr", link1="logit",
                 tau1=seq(.01,.99,.01), condDistobj1=NULL,
                 xformla2=NULL, method2="dr", link2="logit",
                 tau2=seq(.01,.99,.01), condDistobj2=NULL,
                 se=TRUE, iters=100, cl=1) {


    formla <- as.formula(formla)
    dta <- model.frame(terms(formla,data=data),data=data) #or model.matrix
    yname <- colnames(dta)[1]
    tname <- colnames(dta)[2]
    xnames1 <- NULL
    dta1 <- dta
    ##set up the x variables
    if (!(is.null(xformla1))) {
        xformla1 <- as.formula(xformla1)
        xformla1 <- addCovToFormla("0", xformla1)
        xdta1 <- model.matrix(xformla1, data=data)
        dta1 <- cbind.data.frame(dta, xdta1)
        xnames1 <- colnames(xdta1)[-1]
    }
    dta2 <- dta
    xnames2 <- NULL
    ##set up the x variables
    if (!(is.null(xformla2))) {
        xformla2 <- as.formula(xformla2)
        xformla2 <- addCovToFormla("0", xformla2)
        xdta2 <- model.matrix(xformla2, data=data)
        dta2 <- cbind.data.frame(dta, xdta2)
        xnames2 <- colnames(xdta2)[-1]
    }

    cfa1 <- compute.cfa2(tvals, yvals, dta1, yname, tname, xnames1,
                        method1, link1, tau1, condDistobj1)
    cfa2 <- compute.cfa2(tvals, yvals, dta2, yname, tname, xnames2,
                        method2, link2, tau2, condDistobj2)


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
            bdta1 <- dta1[b,]
            bdta2 <- dta2[b,]
            list(bootiter1=compute.cfa2(tvals, yvals, bdta1, yname,
                                       tname, xnames1,
                                       method1, link1, tau1),##$distcondt,
                 bootiter2=compute.cfa2(tvals, yvals, bdta2, yname,
                                       tname, xnames2,
                                       method2, link2, tau2),
                 tvals=bdta1[,tname])
        }, cl=cl)
        bootiterlist1 <- lapply(bstrap, function(x){ x$bootiter1 })
        bootiterlist2 <- lapply(bstrap, function(x){ x$bootiter2 })
        tvallist <- lapply(bstrap, function(x) { x$tvals })
    }

    out <- list(cfa1=CFA.OBJ(tvals, cfa1$distcondt, bootiterlist1, tvallist, coef=cfa1$coef),
                cfa2=CFA.OBJ(tvals, cfa2$distcondt, bootiterlist2, tvallist, coef=cfa2$coef))
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
#' @examples
#' \dontrun{
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' 
#' ## obtain counterfactual results
#' out <- cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr",
#' xformla2=~HEDUC, method2="qr", iters=10, tau1=seq(.05,.95,.05),
#' tau2=seq(.05,.95,.05))
#' test.CFA(out$cfa1, Var, igm$lfincome)
#' }
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


#' @title lige
#'
#' @description compute the local intergenerational elasticity
#'
#' @param cfaobj a CFA object
#' @param h a bandwidth
#' @param se boolean whether or not to compute standard errors
#'
#' @return a CFASE object
#'
#' @examples
#' \dontrun{
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' ## obtain counterfactual results
#' out <- cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr",
#' xformla2=~HEDUC, method2="qr", iters=10, tau1=seq(.05,.95,.05),
#' tau2=seq(.05,.95,.05))
#' lige(out$cfa1, h=0.5)
#' }
#'
#' @export
lige <- function(cfaobj, h, se=T) {
    tvals <- cfaobj$tvals
    bootiterlist <- cfaobj$bootiterlist
    e1res <- getRes.CFA(cfaobj, E, se=FALSE)
    fe1res <- approxfun(tvals, e1res$est)
    lige <- sapply(tvals, function(t) { (fe1res(t+h) - fe1res(t-h))/(2*h) })
    whichT <- which(!is.na(lige))
    lige <- lige[whichT]
    
    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist)) {
            e1resb <- getRes.CFA(bootiterlist[[i]], E, FALSE)
            fe1resb <- approxfun(tvals, e1resb$est)
            ligeb <- sapply(tvals[whichT], function(t) { (fe1resb(t+h) - fe1res(t-h))/(2*h) })
            ##ligeb <- sapply(2:length(theseT), function(i) { (e1resb$est[theseT[i]] - e1resb$est[theseT[i-1]])/(tvals[theseT[i]]-tvals[theseT[i-1]]) } )
            bootout[[i]] <- as.matrix(ligeb)
        }
        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        cb <- unlist(lapply(bootout, function(x) { max(abs((x-lige)/ses)) } ))
        c <- quantile(cb, .95, type=1)
    }
    return(CFASE(tvals=tvals[whichT],est=lige, se=ses, c=c))
}


#' @title Diff.lige
#'
#' @description compute the difference between two estimates of the LIGE
#'
#' @param cfaobj1 the first CFA object
#' @param cfaobj2 the second CFA object
#' @inheritParams lige
#'
#' @return a CFASE object
#' 
#' @examples
#' \dontrun{
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' 
#' ## obtain counterfactual results
#' out <- cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr",
#' xformla2=~HEDUC, method2="qr", iters=10, tau1=seq(.05,.95,.05),
#' tau2=seq(.05,.95,.05))
#' Diff.lige(out$cfa1, out$cfa2, h=0.5)
#' }
#'
#' @export
Diff.lige <- function(cfaobj1, cfaobj2, se=T, h) {
    tvals <- cfaobj1$tvals
    e1res <- getRes.CFA(cfaobj1, E, se=FALSE)
    fe1res <- approxfun(tvals, e1res$est)
    lige1 <- sapply(tvals, function(t) { (fe1res(t+h) - fe1res(t-h))/(2*h) })
    whichT <- which(!is.na(lige1))
    lige1 <- lige1[whichT]

    e2res <- getRes.CFA(cfaobj2, E, se=FALSE)
    fe2res <- approxfun(tvals, e2res$est)
    lige2 <- sapply(tvals, function(t) { (fe2res(t+h) - fe2res(t-h))/(2*h) })
    lige2 <- lige2[whichT]
    
    bootiterlist1 <- cfaobj1$bootiterlist
    bootiterlist2 <- cfaobj2$bootiterlist
    
    out <- lige1 - lige2
    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist1)) {
            e1resb <- getRes.CFA(bootiterlist1[[i]], E, FALSE)
            e2resb <- getRes.CFA(bootiterlist2[[i]], E, FALSE)
            fe1resb <- approxfun(tvals, e1resb$est)
            fe2resb <- approxfun(tvals, e2resb$est)
            lige1b <- sapply(tvals, function(t) { (fe1resb(t+h) - fe1resb(t-h))/(2*h) })
            lige2b <- sapply(tvals, function(t) { (fe2resb(t+h) - fe2resb(t-h))/(2*h) })
            lige1b <- lige1b[whichT]
            lige2b <- lige2b[whichT]
            bootout[[i]] <- as.matrix(lige1b - lige2b)
        }
        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        cb <- unlist(lapply(bootout, function(x) { max(abs((x-out)/ses)) } ))
        c <- quantile(cb, .95, type=1)
    }
    return(CFASE(tvals=tvals[whichT],est=out, se=ses, c=c))
}

#' @title test.lige
#'
#' @description test if the local intergnerational elasticity is the same across
#'  all values of the treatment variable
#'
#' @inheritParams lige
#' @param allt all the values of the treatment variable in the dataset
#'
#' @return a CFASE object
#'
#' @examples
#' \dontrun{
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' 
#' ## obtain counterfactual results
#' out <- cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr",
#' xformla2=~HEDUC, method2="qr", iters=10, tau1=seq(.05,.95,.05),
#' tau2=seq(.05,.95,.05))
#' test.lige(out$cfa1, allt=igm$lfincome, h=0.5)
#' }
#'
#' @export
test.lige <- function(cfaobj, allt, se=T, h) {
    tvals <- cfaobj$tvals
    e1res <- getRes.CFA(cfaobj, E, se=FALSE)
    fe1res <- approxfun(tvals, e1res$est)
    lige <- sapply(tvals, function(t) { (fe1res(t+h) - fe1res(t-h))/(2*h) })
    whichT <- which(!is.na(lige))
    lige <- lige[whichT]
    
    bootiterlist <- cfaobj$bootiterlist##lapply(cfaobj$bootiterlist, function(x) { x$distcondt } )
    tvallist <- cfaobj$tvallist
    
    fcondt <- approxfun(tvals[whichT], lige)
    out <- lige - mean(fcondt(allt), na.rm=TRUE) ## this drops some observations in the extreme tails

    ses <- NULL
    c <- NULL
    if (se) {
        bootout <- list()
        for (i in 1:length(bootiterlist)) {
            e1resb <- getRes.CFA(bootiterlist[[i]], E, FALSE)
            fe1resb <- approxfun(tvals, e1resb$est)
            ligeb <- sapply(tvals, function(t) { (fe1resb(t+h) - fe1resb(t-h))/(2*h) })
            ligeb <- ligeb[whichT]
            bfcondt <- approxfun(tvals[whichT], ligeb)
            bout <- ligeb - mean(bfcondt(tvallist[[i]]), na.rm=T)
            bootout[[i]] <- as.matrix(bout)
        }
        ses <- apply(simplify2array(bootout), c(1,2), sd)
        ses <- if (nrow(ses)==1) as.numeric(ses) else ses
        cb <- unlist(lapply(bootout, function(x) { max(abs((x-out)/ses)) } ))
        c <- quantile(cb, .95, type=1)
    }
    return(CFASE(tvals=tvals[whichT],est=out, se=ses, c=c))
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
#' @examples
#' \dontrun{
#' data(igm)
#' tvals <- seq(10,12,length.out=8)
#' yvals <- seq(quantile(igm$lcfincome, .05), quantile(igm$lcfincome, .95), length.out=50)
#' 
#' ## obtain counterfactual results
#' out <- cfa2(lcfincome ~ lfincome, tvals, yvals, igm, method1="qr",
#' xformla2=~HEDUC, method2="qr", iters=10, tau1=seq(.05,.95,.05),
#' tau2=seq(.05,.95,.05))
#'
#' ## get the difference between the average that adjusts for covariates and
#' ## the one that does not
#' ggplot2.CFA(getResDiff.CFA(out$cfa1, out$cfa2, E), setype="uniform")
#' } 
#'
#' @export
ggplot2.CFA <- function(cfaseobj, setype="pointwise", ylim=NULL,
                        xlabel=NULL, ylabel=NULL, legend=FALSE) {
    tvals <- cfaseobj$tvals
    est <- cfaseobj$est
    se <- cfaseobj$se
    c <- cfaseobj$c
    cmat1 <- cbind.data.frame(tvals, est)
    cmat1$which <- "est"
    if (!is.null(se)) {
        cmat2 <- cbind.data.frame(tvals, se)
        cmat2$which <- "se"
        colnames(cmat2) <- colnames(cmat1)
        cmat <- rbind.data.frame(cmat1, cmat2)
    } else {
        cmat <- cmat1
    }

    cmat <- tidyr::gather(cmat, key=k, value=v, -tvals, -which)

    cmat <- tidyr::spread(cmat, which, v)

    p <- ggplot(cmat, aes(x=tvals, y=est, group=k)) +
        geom_line(aes(color=k))
    if (!is.null(se)) {
        if (setype == "both" | setype=="pointwise") {
            p <- p + geom_line(aes(y=est+1.96*se), lty=3) + geom_line(aes(y=est-1.96*se), lty=3)
        }
        if (setype=="both" | setype=="uniform") {
            p <-  p + geom_line(aes(y=est+c*se), lty=2) + geom_line(aes(y=est-c*se), lty=2)
        }
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

