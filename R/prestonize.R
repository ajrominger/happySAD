#' @title 'Prestonize' an SAD
#'
#' @description Bin an SAD into Preston octaves
#'
#' @details Works for both theoretical SADs (i.e. the density function, by summing
#' probabilty mass) and for empirical SADs (i.e. by binning abundances and summing
#' frequencies)
#'
#' @param x object to Prestonize
#'
#' @return an object of class \code{preston} showing the observed and theoretical
#' probabilities for each Preston octave, and recording the total number of species
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @rdname prestonize
#' @export

prestonize <- function(x) {
    UseMethod('prestonize')
}

#' @rdname prestonize
#' @export

prestonize.numeric  <- function(x) {
    binx <- floor(log(x, 2))
    binx <- table(factor(binx, levels = 0:max(binx)))

    o <- list(binx = as.integer(names(binx)),
              pobs = as.numeric(binx / length(x)),
              pthr = rep(NA, length(binx)),
              S = length(x))

    class(o) <- 'preston'

    return(o)
}

#' @rdname prestonize
#' @export

prestonize.sad <- function(x) {
    if(!is.null(x$data)) {
        obs <- unclass(prestonize(x$data))
    } else {
        obs <- list(S = x$nobs, pobs = NA)
    }

    obs$binx <- 0:13
    obs$pobs <- c(obs$pobs, rep(NA, length(obs$binx) - length(obs$pobs)))


    xx <- 1:(2^max(obs$binx + 1) - 1)
    bb <- floor(log(xx, 2))

    dfun <- getdfun(x)
    pthr <- as.numeric(tapply(dfun(xx), bb, sum))

    o <- list(binx = obs$binx,
              pobs = obs$pobs,
              pthr = pthr,
              S = obs$S)
    class(o) <- 'preston'

    return(o)
}


#' @export
print.preston <- function(x) {
    p <- rbind(x$pobs, x$pthr)
    colnames(p) <- x$binx
    rownames(p) <- c('probs_obs', 'probs_thr')
    print(p)
    cat('S =', x$S, '\n')
}


#' @export
plot.preston <- function(x, ...) {
    sobs <- x$pobs * x$S
    sthr <- x$pthr * x$S

    plot(1, type = 'n',
         xlim = range(x$binx), ylim = c(0, max(sobs, sthr, na.rm = TRUE)),
         xlab = 'Octaves', ylab = 'Number of species', frame.plot = FALSE,
         ...)

    if(any(!is.na(x$pobs))) {
        w <- 0.75 * (0 - par('usr')[1])
        b <- x$binx[!is.na(sobs)]
        sobs <- sobs[!is.na(sobs)]

        rect(xleft = b - w, xright = b + w,
             ybottom = 0, ytop = sobs, col = 'gray')
    }

    points(x$binx, sthr, type = 'b', col = 'red')
}
