dr <- function(y, x, data, yname, xnames) {
    ##form <- as.formula(formla)
    ##dta <- model.frame(terms(formla,data=data),data=data) #or model.matrix
    x <- data.frame(t(x))
    colnames(x) <- xnames
    lgit <- drs.inner(y, data, yname, xnames)
    predict(lgit, newdata=x, type="response")
}

drs <- function(yvals, data, yname, xnames) {
    lapply(yvals, drs.inner, data=data, yname=yname, xnames=xnames)
}

drs.inner <- function(y, data, yname, xnames) {
    IY <- 1*(data[,yname] <= y)
    X <- data[,xnames]
    dta <- cbind.data.frame(IY, X)
    colnames(dta) <- c("IY", xnames)
    formla <- as.formula(paste0("IY ~", paste(xnames, collapse="+")))
    lgit <- glm(formla, data=dta, family=binomial(link=logit))
    lgit
}

distreg <- function(yvals, data, yname, xnames) {
    DR(yvals, drs(yvals, data, yname, xnames))
}

## DR class
DR<- function(yvals, glmlist) {
    out <- list(yvals=yvals, glmlist=glmlist)
    class(out) <- "DR"
    out
}

predict.DR <- function(y, drobj, xdf) {
    yvals <- drobj$yvals
    glmlist <- drobj$glmlist
    if (! (y %in% yvals)) {
         stop("must provide value of y in drobj$yvals")
    }
    x <- xdf
    colnames(x) <-  names(coef(glmlist[[1]]))[-1]
    i <- which(yvals==y)[1]
    predict(glmlist[[i]], newdata=x, type="response")
}
