#' Fit a model with principal components lasso
#'
#' Fit a model using the principal components lasso for an entire regularization
#' path indexed by the parameter \code{lambda}. Fits linear and logistic regression
#' models.
#'
#' The objective function for \code{"gaussian"} is
#'   \deqn{1/2 RSS/nobs + \lambda*||\beta||_1 + \theta/2 \sum quadratic
#'   penalty for group k,}
#' where the sum is over the feature groups \eqn{1, ..., K}. The objective function
#' for \code{"binomial"} is
#'   \deqn{-loglik/nobs + \lambda*||\beta||_1 + \theta/2 \sum quadratic
#'   penalty for group k.}
#'
#' \code{pcLasso} can handle overlapping groups. In this case, the original
#' \code{x} matrix is expanded to a \code{nobs x p_1+...+p_K} matrix (where
#' \code{p_k} is the number of features in group k) such that columns
#' \code{p_1+...+p_{k-1}+1} to \code{p_1+...+p_k} represent the feature matrix for
#' group k. \code{pcLasso} returns the model coefficients for both the expanded
#' feature space and the original feature space.
#'
#' One needs to specify the strength of the quadratic penalty either by
#' specifying \code{ratio}, which is the ratio of shrinkage between the second
#' and first principal components in the absence of the \eqn{\ell_1} penalty,
#' or by specifying the multiplier \code{theta}. \code{ratio} is unitless and is
#' more convenient.
#'
#' \code{pcLasso} always mean centers the columns of the \code{x} matrix. If
#' \code{standardize=TRUE}, \code{pcLasso} will also scale the columns to have
#' standard deviation 1. In all cases, the \code{beta} coefficients returned are
#' for the original \code{x} values (i.e. uncentered and unscaled).
#'
#' @param x Input matrix, of dimension \code{nobs} x \code{nvars}; each row is
#' an observation vector.
#' @param y Response variable. Quantitative for \code{family = "gaussian"}. For
#' \code{family="binomial"}, should be a numeric vector consisting of 0s and 1s.
#' @param w Observation weights. Default is 1 for each observation.
#' @param family Response type. Either \code{"gaussian"} (default) for linear
#' regression or \code{"binomial"} for logistic regression.
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
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of the
#' largest \code{lambda} value. The default depends on the sample size \code{nobs}
#' relative to the number of variables \code{nvars}. If \code{nobs >= nvars},
#' the default is \code{0.0001}, close to zero. If \code{nobs < nvars}, the default is
#' \code{0.01}. This is only used when the user does not specify a \code{lambda}
#' sequence.
#' @param nlam Number of \code{lambda} values; default is 100.
#' @param lambda A user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlam}
#' and \code{lambda.min.ratio}; supplying a value of lambda overrides this.
#' @param standardize If \code{TRUE}, the columns of the feature matrices are
#' standardized before the algorithm is run. Default is \code{FALSE}.
#' @param SVD_info A list containing SVD information. Usually this should not
#' be specified by the user: the function will compute it on its own by default.
#' \code{SVD_info} is a list with three elements:
#' \enumerate{
#'   \item \code{aa}: A row-wise concatenation of the k quadratic penalty matrices.
#'   \item \code{d}: A list where the kth element is the vector of squared singular
#' values of the input matrix corresponding to group k.
#'   \item \code{dd}: A list where the kth element is the vector of
#'   \eqn{d_{k1}^2 - d_{kj}^2} for the input matrix corresponding to group k.
#' }
#' Since the initial SVD of \code{x} can be the largest part of the overall 
#' computation, we allow the user to compute it once and then re-use. See 
#' example below.
#' @param nv Number of singular vectors to use in the singular value
#' decompositions. If not specified, the full SVD is used.
#' @param propack If \code{TRUE} (default), uses \code{propack.svd} from the
#' \code{svd} package to perform the singular value decompositions. If not, uses
#' \code{svd} from base R.
#' @param thr Convergence threhold for the coordinate descent algorithm. Default
#' is \code{1e-4}.
#' @param maxit Maximum number of passes over the data for all lambda values;
#' default is \code{1e5}.
#' @param verbose Print out progess along the way? Default is \code{FALSE}.

#' @return An object of class \code{"pcLasso"}.
#' \item{beta}{If the groups overlap, a \code{p_1+..._p_K x length(lambda)}
#'   matrix of coefficients in the expanded feature space. If not, a \code{nvars}
#'   x \code{length(lambda)} matrix of coefficients in the original feature space.}
#' \item{origbeta}{If the groups overlap, a \code{nvars} x
#'   \code{length(lambda)} matrix of coefficients in the original feature space.
#'   \code{NULL} if the groups do not overlap.}
#' \item{a0}{Intercept sequence of length \code{length(lambda)}.}
#' \item{lambda}{The actual sequence of \code{lambda} values used.}
#' \item{nzero}{If the groups overlap, the number of non-zero coefficients in the
#' expanded feature space for each value of \code{lambda}. Otherwise, the number
#' of non-zero coefficients in the original feature space.}
#' \item{orignzero}{If the groups are overlapping, this is the number of
#'   non-zero coefficients in the original feature space of the model for each
#'   \code{lambda}. If groups are not overlapping, it is \code{NULL}.}
#' \item{jerr}{Error flag for warnings and errors (largely for internal
#'   debugging).}
#' \item{theta}{Value of \code{theta} used in model fitting.}
#' \item{origgroups}{If the \code{groups} parameter was passed to the
#'   function call, this is a copy of that parameter. Otherwise, it is a list of
#'   length 1, with the first element being a vector of integers from 1 to
#'   \code{nvars}.}
#' \item{groups}{If the groups are not overlapping, this has the same
#'   value as \code{groups}. If the groups are overlapping, then
#'   \code{groups[[k]]} is the vector from \code{p_1 + ... + p_{k-1} + 1} to
#'   \code{p_1 + ... p_k}, where \code{p_k} is the number of features in group k.}
#' \item{SVD_info}{A list containing SVD information. See param \code{SVD_info}
#' for more information.}
#' \item{mx}{If groups overlap, column means of the expanded \code{x}
#'   matrix. Otherwise, column means of the original \code{x} matrix.}
#' \item{origmx}{Column means of the original \code{x} matrix.}
#' \item{my}{If \code{family = "gaussian"}, mean of the responses \code{y}.
#'   Otherwise, it is \code{NA}.}
#' \item{overlap}{A logical flag indicating if the feature groups were
#'   overlapping or not.}
#' \item{nlp}{Actual number of passes over the data for all lambda values.}
#' \item{family}{Response type.}
#' \item{call}{The call that produced this object.}
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#'
#' # all features in one group by default
#' fit1 <- pcLasso(x, y, ratio = 0.8)
#' # print(fit1)  # Not run
#' # features in groups
#' groups <- vector("list", 4)
#' for (k in 1:4) {
#'     groups[[k]] <- 5 * (k-1) + 1:5
#' }
#' fit2 <- pcLasso(x, y, groups = groups, ratio = 0.8)
#' # groups can be overlapping
#' groups[[1]] <- 1:8
#' fit3 <- pcLasso(x, y, groups = groups, ratio = 0.8)
#'
#' # specify ratio or theta, but not both
#' fit4 <- pcLasso(x, y, groups = groups, theta = 10)
#'
#' # family = "binomial"
#' y2 <- sample(0:1, 100, replace = TRUE)
#' fit5 <- pcLasso(x, y2, ratio = 0.8, family = "binomial")
#'
#' # example where SVD is computed once, then re-used
#' fit1 <- pcLasso(x, y, ratio = 0.8)
#' fit2 <- pcLasso(x, y, ratio = 0.8, SVD_info = fit1$SVD_info)
#'
#' @export
#' @importFrom svd propack.svd
pcLasso <- function(x, y, w = rep(1,length(y)), family = c("gaussian", "binomial"),
                    ratio = NULL, theta = NULL, groups = vector("list", 1),
                    lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
                    nlam = 100, lambda = NULL, standardize = F, SVD_info = NULL,
                    nv = NULL, propack = T, thr = 1e-4, maxit = 1e5, verbose=FALSE) {
    this.call <- match.call()

    n <- nrow(x)
    p <- ncol(x)
    y <- as.vector(y)
    if (length(y) != n) {
        stop("length of y is not equal to number of rows of x")
    }

    family <- match.arg(family)
    if (family == "binomial" && any(!(unique(y) %in% c(0, 1)))) {
        stop("if family is binomial, y can only contain 0s and 1s")
    }

    # check that each feature appears in some group
    if (length(groups) == 1) groups[[1]] <- 1:p
    if (length(unique(unlist(groups))) < p) {
        stop("Some features not assigned to a group")
    }
    ngroups <- length(groups)
    sizes <- unlist(lapply(groups, length))
    overlap <- F
    origx <- x
    origmx <- colMeans(x)
    origp <- ncol(x)
    origgroups <- groups

    # if there is overlap, redefine x such that features which belong to multiple
    # groups now appear as duplicated columns
    # redefine groups as well (group[[k]] = (p_1+...+p_{k-1}+1):(p_1+...+p_k))
    if (length(origgroups) > 1) {
        nc <- length(unlist(origgroups))
        if (nc > p) {   # there is overlap
            groups <- vector("list", length(origgroups))
            overlap <- T
            x <- matrix(NA, n, nc)
            i1 <- 1
            for(k in 1:ngroups) {
                i2 <- i1 + length(origgroups[[k]]) - 1
                x[, i1:i2] <- origx[, origgroups[[k]]]
                groups[[k]] <- i1:i2
                i1 <- i2 + 1
            }
        }
        p <- ncol(x)
    }

    # scale x and y
    mx <- colMeans(x)
    if (standardize) {
        sx <- apply(x, 2, sd) * sqrt((n-1) / n)
    } else {
        sx <- rep(1, ncol(x))
    }
    x <- scale(x, mx, sx)

    my <- NA
    if (family == "gaussian") {
        my <- mean(y)
        y <- y - my
    }

    # if SVD_info not provided, do the SVD for the feature matrices
    if (is.null(SVD_info)) {
        v = d = dd = vector("list", ngroups)
      if(verbose) cat("Starting SVD computation", fill = T) 
        for (k in 1:ngroups) {
            nvv <- nv
            if (is.null(nv)) {
                nvv <- min(nrow(x[, groups[[k]], drop = F]),
                           ncol(x[, groups[[k]], drop = F]))
            }
            if (nrow(x) > length(groups[[k]])) {
                xtemp <- t(x[, groups[[k]]]) %*% x[, groups[[k]]]
                eig <- eigen(xtemp)
                v[[k]] <- eig$vec
                d[[k]] <- eig$val
                dd[[k]] <- max(d[[k]]) - d[[k]] + .01
            }
            else {
                if (propack) {
                    sv <- svd::propack.svd(x[, groups[[k]], drop = F],
                                           neig = nvv)
                } else {
                    sv <- svd(x[, groups[[k]], drop = F], nv = nvv)
                }
                v[[k]] <- sv$v
                d[[k]] <- sv$d^2
                dd[[k]] <-  max(d[[k]]) - d[[k]] + .01
            }
        }
        if(verbose) cat("SVD completed", fill = T)  
    } else {
        aa <- SVD_info$aa
        d  <- SVD_info$d
        dd <- SVD_info$dd
    }

    # check that exactly one of ratio, theta is provided (also check bounds)
    # if theta is NULL, derive it from ratio
    if (missing(ratio) && missing(theta)) {
        stop("Provide ratio or theta")
    } else if (!missing(ratio) && !missing(theta) &&
               !is.null(ratio) && !is.null(theta)) {
        stop("Provide only ratio or theta, not both")
    } else if (is.null(theta)) {
        if (ratio < 0 || ratio > 1) {
            stop("ratio must be in [0, 1]")
        }
        thetanew <- rep(NA, ngroups)
        for (k in 1:ngroups) {
            thetanew[k] <- d[[k]][2] * (1 - ratio) / (ratio * (d[[k]][1] - d[[k]][2]))
        }
        theta <- mean(thetanew)
    } else if (theta < 0) {
        stop("theta must be non-negative")
    }

    # if SVD_info not provided, compute the VD_{d1^2 - dj^2}V^T matrices for each group
    # and concatenate row-wise
    if (is.null(SVD_info)) {
        aa <- matrix(0, p, max(sizes))
        i1 <- 1
        for(k in 1:ngroups) {
            i2 <- i1 + sizes[k] - 1
            aa[i1:i2, 1:sizes[k]] <- scale(v[[k]], FALSE, 1/(dd[[k]])) %*% t(v[[k]])
            i1 <- i2 + 1
        }
        # build up SVD_info for returned output
        SVD_info <- list()
        SVD_info$aa <- aa
        SVD_info$d  <- d
        SVD_info$dd <- dd
    }

    # compute mg: vector of length length(sizes)+1. Denotes the number for the
    # first feature in group k, with an additional last entry = p + 1.
    i1 <- 1
    mg <- c(1, cumsum(sizes) + 1)

    # get lambda sequence
    ulam <- lambda * n   # we need this multiplication by n because documentation has RSS/(2n),
                         # but fortran code assumes RSS/2
    if (is.null(lambda)) {
        maxlam <- max(abs(t(x) %*% (y - mean(y))))
        if (family == "binomial") {
            maxlam <- 4 * maxlam
        }
        ulam <- exp(seq(log(maxlam), log(maxlam * lambda.min.ratio), length = nlam))
    }
    nlam <- length(ulam)

    out <- pcLassof(x, y, w, mg, aa, ulam, theta, ngroups, thr=thr, maxit=maxit,
                    family=family, verbose=verbose)
    
    nzero <- colSums(out$beta != 0)
    
    # get intercept sequence
    if (family == "gaussian") {
        a0 <- rep(my, nlam)
    } else {
        # if here, family must be binomial
        a0 <- out$a0
    }
    
    # rescale intercept values
    a0 <- a0 - colSums(out$beta * mx / sx)

    # rescale beta values (beta will only change when standardize = TRUE)
    out$beta <- t(scale(t(out$beta), F, sx))

    # convert beta values in expanded space to beta values in original space
    origbeta <- NULL
    orignzero <- NULL
    if (overlap) {
        origbeta <- matrix(0, ncol(origx), ncol(out$beta))
        for(k in 1:ngroups) {
            origbeta[origgroups[[k]], ] <- origbeta[origgroups[[k]], ] +
                out$beta[groups[[k]], ]
        }
        orignzero <- colSums(origbeta != 0)
    }

    lambda <- out$ulam / n
    if (family == "binomial") {
        lambda=lambda / 4
    }
    out <- list(beta=out$beta, origbeta=origbeta, a0=a0, lambda=lambda,
                nzero=nzero, orignzero=orignzero, jerr=out$jerr, theta=theta,
                origgroups=origgroups, groups=groups, SVD_info = SVD_info, mx=mx,
                origmx=origmx, my=my, overlap=overlap, nlp=out$nlp,
                family=family, call=this.call)
    
    yhat <- predict.pcLasso(out, origx)
    if (family == "gaussian") dev <- colSums(apply(yhat,2,msefun, y))
    if (family == "binomial") dev <- colSums(apply(yhat,2,binfun, y))
    out$dev <- dev

    class(out) <- "pcLasso"
    return(out)
}
