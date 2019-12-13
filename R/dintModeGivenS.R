#' @title Mode probabilities
#'
#' @description Probability of a mode across Preston octaves
#'
#' @details For Preston octaves 0 through 12, this function calculates the probability that each
#' octave could be a modal octave (mode defined as equal to the maximum frequency of species).
#' More than one octave could meet this criterian, so the returned vector of probabilities will
#' not neccesarily sum to 1.  The probability of a mode at octave i is
#' \eqn{\sum_s P(o_i = s) P(\text{all } o_{j \neq i} \leq s)} where \eqn{P(o_i = s)} is simply
#' \code{dbinom(s, S, p[i])} where \code{S} is the total number of species and \code{p[i]} is the
#' summed probability of the SAD corresponding to octave i. And where
#' \eqn{P(\text{all } o_{j \neq i} \leq s)} is given by \code{pmultinom(rep(s, 12), S - s, p[-i])}
#'
#' @param S number of spcies
#' @param mod function specifying the SAD density function (e.g. \code{dfish})
#' @param pars vector of parameters for the SAD model
#'
#' @return a vector of probabilities equal in length to number of octaves (13)
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

dintModeGivenS <- function(S, mod, pars) {
    # browser()
    noct <- 13
    x <- 1:noct
    nn <- 1:(2^(noct) - 1)
    octs <- floor(log(nn, 2))
    pp <- do.call(mod, c(list(nn), as.list(pars)))
    octProbs <- tapply(pp, octs, sum)


    possModes <- floor(S / noct):S

    # limit number of octaves we look at
    goodx <- sapply(octProbs, function(p) {
        dbinom(possModes, S, p) > .Machine$double.eps^0.75
    })
    goodx <- colSums(goodx) > 0

    ii <- range(which(dbinom(possModes, S, max(octProbs[goodx])) > .Machine$double.eps^0.75),
                which(dbinom(possModes, S, min(octProbs[goodx])) > .Machine$double.eps^0.75))

    # limit the number of modes we're integrating over
    likModes <- possModes[ii[1]:ii[2]]


    # probability of S in bin 1
    probSinBi <- dbinom(likModes, S, octProbs[1])

    # prob all other bins have less
    probRestLessS <- numeric(length(probSinBi))

    # only calculate costly convolution for cases where it matters
    probRestLessS[probSinBi > .Machine$double.eps] <-
        sapply(likModes[probSinBi > .Machine$double.eps], function(thisS) {
            o <- pmultinom(rep(thisS, noct - 1), S - thisS, octProbs[-1])
            o[o > 1] <- 1
            o[o < .Machine$double.eps] <- 0

            return(o)
        })

    # sum over all S
    allS <- probSinBi * probRestLessS

    # return prob that first bin is not modal
    return(1 - sum(allS))


    #
    #
    # # browser()
    # o <- sapply(x[goodx], function(i) {
    #     probSinBi <- dbinom(likModes, S, octProbs[i])
    #     probRestLessS <- numeric(length(probSinBi))
    #
    #     # only calculate costly convolution for cases where it matters
    #     probRestLessS[probSinBi > .Machine$double.eps] <-
    #         sapply(likModes[probSinBi > .Machine$double.eps], function(thisS) {
    #             o <- pmultinom(rep(thisS, noct - 1), S - thisS, octProbs[-i])
    #             o[o > 1] <- 1
    #             o[o < .Machine$double.eps] <- 0
    #
    #             return(o)
    #         })
    #
    #     allS <- probSinBi * probRestLessS
    #     return(sum(allS))
    # })
    #
    # out <- numeric(noct)
    # out[goodx] <- o
    #
    # return(out)
}

appx_dintModeGivenS <- function(S, pars, mod) {
    noct <- 13
    x <- 1:noct
    nn <- 1:(2^(noct) - 1)
    octs <- floor(log(nn, 2))
    pp <- do.call(mod, c(list(nn), as.list(pars)))
    octProbs <- tapply(pp, octs, sum)


    possModes <- floor(S / noct):S

    # limit number of octaves we look at
    goodx <- sapply(octProbs, function(p) {
        dbinom(possModes, S, p) > .Machine$double.eps^0.75
    })
    goodx <- colSums(goodx) > 0

    ii <- range(which(dbinom(possModes, S, max(octProbs[goodx])) > .Machine$double.eps^0.75),
                which(dbinom(possModes, S, min(octProbs[goodx])) > .Machine$double.eps^0.75))

    # limit the number of modes we're integrating over
    likModes <- possModes[ii[1]:ii[2]]

    # browser()
    o <- sapply(x[goodx], function(i) {
        probSinBi <- dbinom(likModes, S, octProbs[i])
        probRestLessS <- numeric(length(probSinBi))

        # only calculate costly convolution for cases where it matters
        probRestLessS[probSinBi > .Machine$double.eps] <-
            sapply(likModes[probSinBi > .Machine$double.eps], function(thisS) {
                # browser()
                pmultinomapprox(rep(thisS, noct - 1), S - thisS, octProbs[-i])
            })

        allS <- probSinBi * probRestLessS
        return(sum(allS))
    })

    out <- numeric(noct)
    out[goodx] <- o

    return(out)
}
