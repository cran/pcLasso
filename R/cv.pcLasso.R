#' Cross-validation for pcLasso
#'
#' Does \code{k}-fold cross-validation for \code{pcLasso}.
#'
#' This function runs \code{pcLasso nfolds+1} times: the first to get the
#' \code{lambda} sequence, and the remaining \code{nfolds} times to compute the
#' fit with each of the folds omitted. The error is accumulated, and the mean
#' error and standard deviation over the folds is compued. Note that
#' \code{cv.pcLasso} does NOT search for values of \code{theta} or \code{ratio}.
#' A specific value of \code{theta} or \code{ratio} should be supplied.
#'
#' @param x \code{x} matrix as in \code{pcLasso}.
#' @param y \code{y} matrix as in \code{pcLasso}.
#' @param w Observation weights. Default is 1 for each observation.
#' @param ratio Ratio of shrinkage between the second and first principal components
#' in the absence of the \eqn{\ell_1} penalty. More convenient way to specify the
#' strength of the quadratic penalty. A value between 0 and 1 (only 1 included).
#' \code{ratio = 1} corresponds to the lasso, 0.5-0.9 are good values to use.
#' Default is \code{NULL}. Exactly one of \code{ratio} or \code{theta} must be
#' specified.
#' @param theta Multiplier for the quadratic penalty: a non-negative real number.
#' \code{theta = 0} corresponds to the lasso, and larger \code{theta} gives
#' strong shrinkage toward the top principal components. Default is \code{NULL}.
#' Exactly one of \code{ratio} or \code{theta} must be specified.
#' @param groups A list describing which features belong in each group. The
#' length of the list should be equal to the number of groups, with
#' \code{groups[[k]]} being a vector of feature/column numbers which belong to
#' group k. Each feature must be assigned to at least one group. Features can
#' belong to more than one group. By default, all the features belong to a
#' single group.
#' @param family Response type. Either \code{"gaussian"} (default) for linear
#' regression or \code{"binomial"} for logistic regression.
#' @param nfolds Number of folds for CV (default is 10). Although \code{nfolds}
#' can be as large as the sample size (leave-one-out CV), it is not recommended
#' for large datasets. Smallest value allowable is \code{nfolds = 3}.
#' @param foldid An optional vector of values between 1 and \code{nfold}
#' identifying what fold each observation is in. If supplied, \code{nfold} can
#' be missing.
#' @param keep If \code{keep = TRUE}, a prevalidated array is returned
#' containing fitted values for each observation at each value of lambda. This
#' means these fits are computed with this observation and the rest of its fold
#' omitted. Default is \code{FALSE}.
#' @param verbose Print out progess along the way? Default is \code{FALSE}.
#' @param ... Other arguments that can be passed to \code{pcLasso}.
#'
#' @return An object of class \code{"cv.pcLasso"}, which is a list with the
#' ingredients of the cross-validation fit.
#' \item{glmfit}{A fitted \code{pcLasso} object for the full data.}
#' \item{theta}{Value of \code{theta} used in model fitting.}
#' \item{lambda}{The values of \code{lambda} used in the fits.}
#' \item{nzero}{If the groups overlap, the number of non-zero coefficients
#'   in the model \code{glmfit} for the expanded feature space, at each value of
#'   \code{lambda}. Otherwise, the number of non-zero coefficients in the model
#'   \code{glmfit} for the original feature space.}
#' \item{orignzero}{If the groups are overlapping, this is the number of
#'   non-zero coefficients in the model \code{glmfit} for the original feature
#'   space, at each \code{lambda}. If groups are not overlapping, it is
#'   \code{NULL}.}
#' \item{fit.preval}{If \code{keep=TRUE}, this is the array of prevalidated
#'   fits.}
#' \item{cvm}{The mean cross-validated error: a vector of length
#'   \code{length(lambda)}.}
#' \item{cvse}{Estimate of standard error of \code{cvm}.}
#' \item{cvlo}{Lower curve = \code{cvm - cvsd}.}
#' \item{cvup}{Upper curve = \code{cvm + cvsd}.}
#' \item{lambda.min}{The value of \code{lambda} that gives minimum
#'   \code{cvm}.}
#' \item{lambda.1se}{The largest value of \code{lambda} such that the CV
#'   error is within one standard error of the minimum.}
#' \item{foldid}{If \code{keep=TRUE}, the fold assignments used.}
#' \item{name}{Name of error measurement used for CV.}
#' \item{call}{The call that produced this object.}
#'
#' @seealso \code{\link{pcLasso}} and \code{\link{plot.cv.pcLasso}}.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' groups <- vector("list", 4)
#' for (k in 1:4) {
#'     groups[[k]] <- 5 * (k-1) + 1:5
#' }
#' cvfit1 <- cv.pcLasso(x, y, groups = groups, ratio = 0.8)
#'
#' # change no. of CV folds
#' cvfit2 <- cv.pcLasso(x, y, groups = groups, ratio = 0.8, nfolds = 5)
#' # specify which observations are in each fold
#' foldid <- sample(rep(seq(5), length = length(y)))
#' cvfit3 <- cv.pcLasso(x, y, groups = groups, ratio = 0.8, foldid = foldid)
#'
#' # keep=TRUE to have pre-validated fits and foldid returned
#' cvfit4 <- cv.pcLasso(x, y, groups = groups, ratio = 0.8, keep = TRUE)
#'
#' @export
cv.pcLasso <- function(x, y, w = rep(1,length(y)), ratio = NULL, theta = NULL,
                       groups = vector("list", 1), family = "gaussian",
                       nfolds = 10, foldid = NULL, keep = FALSE, verbose=FALSE, ...) {
    this.call <- match.call()
    n <- nrow(x)
    p <- ncol(x)

    if(length(groups) == 1) groups[[1]] <- 1:p
    ngroups <- length(groups)

    # get CV fold information
    if (!missing(foldid)) {
        nfolds <- length(unique(foldid))
    } else {
        foldid <- sample(rep(seq(nfolds), length = length(y)))
    }
    if (nfolds < 3) {
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    }

    # get pcLasso fit for all of the data
    fit0 <- pcLasso(x, y, groups=groups, ratio=ratio, theta=theta,
                    family=family, ...)
    cat("Initial fit done- including SVD", fill = T)
    fits <- vector("list", nfolds)
    for (ii in 1:nfolds) {
        cat(c("Fold=",ii), fill = T)
        oo <- foldid == ii
        xc <- x[!oo, , drop = F]
        yy <- y[!oo]
        fits[[ii]] <- pcLasso(xc, yy, SVD_info = fit0$SVD_info,
                              groups=groups, theta=fit0$theta,
                              family=family, lambda=fit0$lambda, ...)
    }

    # get predictions
    yhat <- matrix(NA, n, length(fit0$lambda))
    for (ii in 1:nfolds) {
        oo <- foldid == ii
        out <- predict(fits[[ii]], x[oo, , drop = F])
        yhat[oo, 1:ncol(out)] <- out
    }
    if (family == "binomial") {
        yhat <- 1 / (1 + exp(-yhat))
    }

    if (family == "gaussian") {
        errfun = msefun
        name = "Mean-Squared Error"
    }
    if (family == "binomial") {
        errfun = binfun
        name = "Deviance"
    }

    # compute CV error
    ym <- array(y, dim(yhat))
    err <- errfun(yhat,ym)
    cvm <- apply(err, 2, mean, na.rm = T)
    nn <- apply(!is.na(err), 2, sum, na.rm = T)
    cvse <- sqrt(apply(err, 2, var, na.rm = T) / nn)
    cvlo <- cvm - cvse
    cvup <- cvm + cvse

    yhat.preval <- NULL
    foldid_copy <- NULL
    if (keep) {
        yhat.preval <- yhat
        foldid_copy <- foldid
    }

    # get best and 1se lambda
    imin <- which.min(cvm)
    lambda.min <- fit0$lambda[imin]
    imin.1se <- which(cvm < cvm[imin] + cvse[imin])[1]
    lambda.1se <- fit0$lambda[imin.1se]

    obj <- list(glmfit=fit0, theta=fit0$theta, lambda=fit0$lambda,
                nzero=fit0$nzero, orignzero=fit0$orignzero, fit.preval=yhat,
                cvm=cvm, cvse=cvse, cvlo=cvlo, cvup=cvup, lambda.min=lambda.min,
                lambda.1se=lambda.1se, foldid=foldid_copy, name=name,
                call=this.call)
    class(obj) <- "cv.pcLasso"
    return(obj)
}
