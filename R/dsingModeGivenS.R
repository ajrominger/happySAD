#' @title Probability of a single mode
#'
#' @description For a given model and number of species, returns the probability of a
#' single mode in the Preston plot
#'
#' @details A mode is defined by the cuttoff \code{modeCutOff}. This can be made a proportional cuttoff if \code{relative = TRUE} is relative
#' to the height of the peak. For example if \code{modeCutOff = 0.1, relative = TRUE} this sequence would contain a mode:
#' \code{c(17, 20, 17)} but this one would not: \code{c(19, 20, 19)}.  The Probability of a single mode
#' is numerically involved and so a Monte Carlo approximation is used with \code{B} replicates.
#'
#' @param S number of species
#' @param mod the SAD model specified by the \code{d*} function, e.g. \code{dfish}
#' @param pars model parameters
#' @param modeCutOff the fraction by which a peak must exceed its neigbors to count as a mode
#' @param relative logical, should \code{modeCutOff} be taken as relative or not
#' @param B number of Monte Carlo iterations
#'
#' @return a numeric giving the desired probability
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

dsingModeGivenS <- function(S, mod, pars, modeCutOff = 0.1, relative = TRUE, B = 1e+04) {
    noct <- 13
    x <- 1:noct
    nn <- 1:(2^(noct) - 1)
    octs <- floor(log(nn, 2))
    pp <- do.call(mod, c(list(nn), as.list(pars)))
    octProbs <- tapply(pp, octs, sum)

    r <- rmultinom(B, S, octProbs)

    o <- lapply(1:ncol(r), function(i) {
        return(.allLessMax(r[, i], modeCutOff, relative))
        return(.nmode(r[, i], modeCutOff) == 1)
    })

    return(mean(unlist(o)))
}

.allLessMax <- function(x, cutOff, relative) {
    imax <- which.max(rle(x)$values)

    if(relative) cutOff <- max(x) * cutOff

    maxsplit <- list(x[1:imax], x[imax:length(x)])

    return(all(diff(maxsplit[[1]]) > 1 - cutOff) &
               all(diff(maxsplit[[2]]) < -1 + cutOff))
}

#
# .nmode <- function(x, cutOff) {
#     y <- log(x)
#     y[y < 0] <- -1
#
#     length(localMaxima(y, by = log(1 / (1 - cutOff)), xmin = -1))
# }
