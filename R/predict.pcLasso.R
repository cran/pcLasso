#' Make predictions from a "pcLasso" object
#'
#' This function returns the predictions from a "\code{pcLasso}" object
#' for a new data matrix. 
#'
#' Note that \code{xnew} should have the same number of columns as the original
#' feature space, regardless of whether the groups are overlapping or not.
#'
#' @param object Fitted "\code{pcLasso}" object.
#' @param xnew Matrix of new values for \code{x} at which predictions are to
#' be made.
#' @param ... Potentially other arguments to be passed to and from methods;
#' currently not in use.
#'
#' @return Predictions of \eqn{E(y|xnew)} which the model \code{object} makes at 
#' \code{xnew}. These are probabilities for the binomial family.
#'
#' @seealso \code{\link{pcLasso}}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#'
#' # family = "gaussian"
#' y <- rnorm(100)
#' fit1 <- pcLasso(x, y, ratio = 0.8)
#' predict(fit1, xnew = x[1:5, ])
#'
#' # family = "binomial"
#' y2 <- sample(0:1, 100, replace = TRUE)
#' fit2 <- pcLasso(x, y2, ratio = 0.8, family = "binomial")
#' predict(fit2, xnew = x[1:5, ])
#'
#' @export
predict.pcLasso <- function(object, xnew, ...) {
    if (object$overlap) {
        beta <- object$origbeta
    } else {
        beta <- object$beta
    }

    out <- t(object$a0 + t(xnew %*% beta))
    if (object$family == "binomial") {
        out <- 1 / (1 + exp(-out))
    }
    return(out)
}
