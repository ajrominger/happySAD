#' @title Multinomial CDF
#'
#' @description Cumulative probability function of the multinomial distribution
#'
#' @details For Preston octaves 0 through 12, this function calculates the probability that each
#' octave could be a modal octave (mode defined as equal to the maximum frequency of species).
#' More than one octave could meet this criterian, so the returned vector of probabilities will
#' not neccesarily sum to 1.  The probability of a mode at octave i is
#' $\sum_s P(o_i = s) P(\text{all } o_{j \neq i} \leq s)$ where $P(o_i = s)$ is simply
#' \code{dbinom(s, S, p[i])} where \code{S} is the total number of species and \code{p[i]} is the
#' summed probability of the SAD corresponding to octave i. And where
#' $P(\text{all } o_{j \neq i} \leq s)$ is given by \code{pmultinom(rep(s, 12), S - s, p[-i])}
#'
#' @param S number of spcies
#' @param pars vector of parameters for the SAD model
#' @param mod string specifying the SAD density function (e.g. \code{'dfish'})
#'
#' @return a vector of probabilities equal in length to number of octaves (13)
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

dintModeGivenS <- function(S, pars, mod) {
    noct <- 13
    x <- 1:noct
    nn <- 1:(2^(noct) - 1)
    octs <- floor(log(nn, 2))
    pp <- do.call(mod, c(list(nn), as.list(pars)))
    octProbs <- tapply(pp, octs, sum)


    possModes <- ceiling(S / noct):S
    ii <- range(which(dbinom(possModes, S, max(octProbs)) > .Machine$double.eps),
                which(dbinom(possModes, S, min(octProbs)) > .Machine$double.eps))

    # limit the number of things we're integrating over
    likModes <- possModes[ii[1]:ii[2]]
    goodx <- dbinom(min(likModes), S, octProbs) > .Machine$double.eps

    o <- sapply(x[goodx], function(i) {
        allS <- sapply(likModes, function(thisS) {
            probSinBi <- dbinom(thisS, S, octProbs[i])
            probRestLessS <- pmultinom(rep(thisS, noct - 1), S - thisS, octProbs[-i])
            return(probSinBi * probRestLessS)
        })

        return(sum(allS))
    })

    out <- numeric(noct)
    out[goodx] <- o

    return(out)
}

