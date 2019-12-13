#' @title Exponential axis
#'
#' @description When plotting already log-transformed data, add an axis with exponential
#' tick labels
#'
#' @details assumes \eqn{log_{10}}
#'
#' @param side side to add axis to
#' @param at locations of ticks
#' @param ... additional parameters passed to \code{axis}
#'
# @return
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

expAxis <- function(side, at, ...) {
    labs <- sapply(at, function(p) {
        eval(substitute(expression(10^p), list(p = p)))
    })

    axis(side, at, labels = labs, ...)
}
