# helper functions for pcLasso, no need to export
pcLassof <- function(x, y, w, mg, aa, ulam, theta, ngroups, thr = 1e-4,
                     maxit = 1e5, family = "gaussian", verbose = FALSE) {
    no = nrow(x)
    ni = ncol(x)
    ne = ni
    nx = ni
    ng = ngroups
    if (is.null(w)) w = rep(1,no)
    nlam=length(ulam)
    verbose=1*verbose
    mode(no)="integer"
    mode(ni)="integer"
    mode(x)="double"
    mode(y)="double"
    mode(w)="double"
    mode(theta)="double"
    mode(ng)="integer"
    mode(mg)="integer"
    mode(aa)="double"
    mode(ne)="integer"
    mode(nx)="integer"
    mode(nlam)="integer"
    mode(ulam)="double"
    mode(thr)="double"
    mode(maxit)="integer"
  mode(verbose)="integer"
    if (family == "gaussian") {
        out = .Fortran("pclasso",
                       no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam=nlam,ulam=ulam,thr,maxit,verbose,
                       ao=double(nx*nlam),
                       ia=integer(nx),
                       kin=integer(nlam),
                       nlp=integer(1),
                       jerr=integer(1),
                       PACKAGE="pcLasso")
    }
    if (family == "binomial") {
        out = .Fortran("logpclasso",
                       no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam=nlam,ulam=ulam,thr,maxit,verbose,
                       a0=double(nlam),
                       ao=double(nx*nlam),
                       ia=integer(nx),
                       kin=integer(nlam),
                       nlp=integer(1),
                       jerr=integer(1),
                        PACKAGE="pcLasso")
    }

    ao = matrix(out$ao, nrow = ni)

    # uncompress soln
    beta <- matrix(0, ni, out$nlam)
    for (klam in 1:out$nlam) {
        temp <- out$kin[klam]
        beta[out$ia[1:temp], klam] <- ao[1:temp, klam]
    }

    a0 <- NA
    if (family == "binomial") a0 <- out$a0

    return(list(beta=beta, a0=a0, ulam=out$ulam, nlam=nlam, nlp=out$nlp,
                jerr=out$jerr))
}

msefun <- function(yhat,y) {
    (y - yhat)^2
}

binfun <- function(yhat, y) {
    - y * log(yhat) - (1 - y) * log(1 - yhat)
}

error.bars <- function(x, upper, lower, width = 0.02, ...) {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}


print.pcLasso=function(x,digits = max(3, getOption("digits") - 3),...){
       devratio=(x$dev[1]-x$dev)/x$dev[1]
  
cat("\nCall: ", deparse(x$call), "\n\n")
    print(cbind(Nonzero = x$nzero, `%Dev` = signif(devratio, digits), 
        Lambda = signif(x$lambda, digits)))

                
}
