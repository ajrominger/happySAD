#' @title Expected N given S and an SAD model
#'
#' @description The expected value of N, the sum of all abundances over all species
#'
#' @details Uses the \code{dNGivenS} function internally.
#'
#' @param S number of spcies
#' @param N total number of individuals
#' @param mod string specifying the SAD model (e.g. \code{'fish'})
#' @param pars vector of parameters for the SAD model
#'
#' @return \code{expectedN} returns the expectation of N, while \code{solve4expectedN}
#' finds the value of S that gives the desired expected value of N given the SAD model
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export
#' @rdname expectedN

expectedN <- function(S, mod, pars) {
    N <- min(S):(1000 * max(S))
    xp <- sapply(S, function(s) {
        n <- N[N > s]
        o <- dNGivenS(n, s, pars, mod) * n
        if(o[length(o)] > .Machine$double.eps^0.5) {
            cat('need more \n')
            nplus <- (max(n) + 1):(max(n) * 1000)
            oplus <- dNGivenS(nplus, s, pars, mod) * nplus
            o <- c(o, oplus)
        }

        return(o)
    })

    sum(xp)
}


#' @export
#' @rdname expectedN

solve4expectedN <- function(N, mod, pars) {
    uniroot(.solExpectedN, c(1, N), N = N, mod = 'fish', pars = bta)
}


.solExpectedN <- function(S, N, mod, pars) {
    expectedN(S, mod, pars) - N
}
