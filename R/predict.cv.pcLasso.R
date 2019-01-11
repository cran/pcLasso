#' Make predictions from a "cv.pcLasso" object
#'
#' This function returns the predictions for a new data matrix from a
#' cross-validated pcLasso model by using the stored "\code{glmfit}" object and
#' the optimal value chosen for \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction. Note that \code{xnew} should have the same number of columns as
#' the original feature space, regardless of whether the groups are overlapping
#' or not.
#'
#' @param object Fitted "\code{cv.pcLasso}" object.
#' @param xnew Matrix of new values for \code{x} at which predictions are to
#' be made.
#' @param s Value of the penalty parameter \code{lambda} at which predictions are
#' required. Default is the value \code{s="lambda.1se"} stored in the CV
#' \code{fit}. Alternatively, \code{s="lambda.min"} can be used.
#' @param ... Potentially other arguments to be passed to and from methods;
#' currently not in use.
#'
#' @return Predictions which the cross-validated model makes for \code{xnew} at
#' the optimal value of \code{lambda}. Note that the default is the "lambda.1se" for lambda,
#' to make this function consistent with \code{cv.glmnet} in the \code{glmnet} 
#' package. The output is predictions of \eqn{E(y|xnew)}: these are probabilities 
#' for the binomial family.
#'
#' @seealso \code{\link{cv.pcLasso}} and \code{\link{predict.pcLasso}}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#'
#' cvfit <- cv.pcLasso(x, y, ratio = 0.8)
#' predict(cvfit, xnew = x[1:5, ])
#' predict(cvfit, xnew = x[1:5, ], s = "lambda.min")
#'
#' @export
predict.cv.pcLasso <- function(object, xnew, s = c("lambda.1se", "lambda.min"),
                               ...) {
    s <- match.arg(s)
    lambda <- object[[s]]
    idx <- which(object$lambda == lambda)

    predictions <- predict(object$glmfit, xnew)
    return(predictions[, idx])
}
