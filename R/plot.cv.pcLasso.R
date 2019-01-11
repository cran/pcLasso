#' Plot the cross-validation curve produced by "cv.pcLasso" object
#'
#' Plots the cross-validation curve produced by a \code{cv.pcLasso} object, along
#' with upper and lower standard deviation curves, as a function of the \code{lambda}
#' values used.
#'
#' A plot is produced and nothing is returned.
#'
#' @param x Fitted "\code{cv.pcLasso}" object.
#' @param sign.lambda Either plot against \code{log(lambda)} (default) or
#' \code{-log(lambda)} (if \code{sign.lambda = -1}).
#' @param orignz If \code{TRUE} (default), prints the number of non-zero
#' coefficients in the original feature space. If not, prints the number of
#' non-zero coefficients in the expanded feature space. No effect if groups are
#' not overlapping.
#' @param ... Other graphical paramters to plot.
#' @seealso \code{\link{pcLasso}} and \code{\link{cv.pcLasso}}.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' groups <- vector("list", 4)
#' for (k in 1:4) {
#'     groups[[k]] <- 5 * (k-1) + 1:5
#' }
#' cvfit <- cv.pcLasso(x, y, ratio = 0.8, groups = groups)
#' plot(cvfit)
#' # plot flipped: x-axis tracks -log(lambda) instead
#' plot(cvfit, sign.lambda = -1)
#'
#' # if groups overlap, orignz can be used to decide which space to count the
#' # number of non-zero coefficients at the top
#' groups[[1]] <- 1:8
#' cvfit <- cv.pcLasso(x, y, ratio = 0.8, groups = groups)
#' plot(cvfit)                  # no. of non-zero coefficients in original space
#' plot(cvfit, orignz = FALSE)  # no. of non-zero coefficients in expanded space
#'
#' @export
plot.cv.pcLasso <- function(x, sign.lambda = 1, orignz = TRUE, ...) {
    xlab <- "log(Lambda)"
    if (sign.lambda < 0) {
        xlab = paste("-", xlab, sep = "")
    }
    plot.args <- list(x = sign.lambda * log(x$lambda), y = x$cvm,
                      ylim = range(x$cvup, x$cvlo),
                      xlab = xlab, ylab = x$name, type = "n")
    new.args <- list(...)
    if (length(new.args)) {
        plot.args[names(new.args)] = new.args
    }
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(x$lambda), x$cvup, x$cvlo,
               width = 0.01, col = "darkgrey")
    points(sign.lambda * log(x$lambda), x$cvm, pch = 20, col = "red")

    # what to print at top depends on whether groups are overlapping and orignz value
    if (is.null(x$orignzero) || !orignz) {  # groups are non-overlapping
        top_labels <- x$nzero
    } else {
        top_labels <- x$orignzero
    }
    axis(side = 3, at = sign.lambda * log(x$lambda), labels = paste(top_labels),
         tick = FALSE, line = 0)

    abline(v = sign.lambda * log(x$lambda.min), lty = 3)
    abline(v = sign.lambda * log(x$lambda.1se), lty = 3)
    invisible()
}
