#' @title Species by beta image
#'
#' @description An image plot customized for a matrix of probabilities across different
#' numbers of species and beta parameters
#'
#' @details assumes \code{b} was made over a \eqn{log_{10}} scale
#'
#' @param s the species values used to calculate \code{p}
#' @param b the beta values used to calculate \code{p}
#' @param p the matrix of probabilities (rows correspond to species values)
#' @param colz vector of colors to use
#' @param list of arguments to pass to \code{par}, note, at a minimum must include \code{mar}
#' @param legendLab lable for the legend
#' @param ... additional parameters passed to \code{image}
#'
# @return
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

sbimage <- function(s, b, p, colz, parArgs, legendLab, ...) {
    par(parArgs)

    layout(matrix(1:2, nrow = 1), widths = c(3, 0.75))
    image(log(b, 10), s, t(p), col = colz, zlim = c(0, 1),
          xlab = expression(beta), ylab = 'S', xaxt = 'n')
    expAxis(1, pretty(log(b, 10)))

    par(mar = c(parArgs$mar[1], parArgs$mar[4], parArgs$mar[3], parArgs$mar[2]))
    plot(1, ylim = c(0, 1), axes = FALSE, type = 'n', xlab = '', ylab = '')

    rect(xleft = par('usr')[1], xright = par('usr')[2],
         ybottom = seq(0, 1, length.out = length(colz) + 1)[- (length(colz) + 1)],
         ytop = seq(0, 1, length.out = length(colz) + 1)[-1],
         col = colz, border = colz)
    rect(xleft = par('usr')[1], xright = par('usr')[2], ybottom = 0, ytop = 1)
    axis(4)
    mtext(legendLab, side = 4, line = 1.75)
}


