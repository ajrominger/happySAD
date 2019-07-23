#' @title Probability of N given S
#'
#' @description Density, distribution, and quantile functions for the total number of
#' individuals given the number of species
#'
#' @details These functions assume a sufficiently large number of species (in practice >= 20)
#' such that the Central Limit Theorem applies and we can approxiate the distribution of the
#' sum of \code{S} random variables as a normal distribution
#'
#' @param N total abundance of interest
#' @param S number of species
#' @param p probability for which to return the quantile
#' @param pars parameter(s) for the SAD model
#' @param mod a string specifying the SAD model
#' @param log logical, should the log probability be used
#' @param log.p logical, should the log cumulative probability be used
#' @param lower.tail logical, should the lower tail be used
#'
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export
#' @rdname NGivenS

dNGivenS <- function(N, S, pars, mod, log = FALSE) {
    musig <- switch(mod,
                    'fish' = list(.fishMean(pars) * S, sqrt(.fishVar(pars) * S)),
                    'plnorm' = list(0, 0),
                    'stick' = list(0, 0),
                    'tnegb' = list(0, 0),
                    'untb' = list(0, 0))

    return(dnorm(N, musig[[1]], musig[[2]], log = log))
}

#' @export
#' @rdname NGivenS
pNGivenS <- function(N, S, pars, mod, lower.tail = TRUE, log.p = FALSE) {
    musig <- switch(mod,
                    'fish' = list(.fishMean(pars) * S, sqrt(.fishVar(pars) * S)),
                    'plnorm' = list(0, 0),
                    'stick' = list(0, 0),
                    'tnegb' = list(0, 0),
                    'untb' = list(0, 0))

    return(pnorm(N, musig[1], musig[2], lower.tail = lower.tail, log.p = log.p))
}

#' @export
#' @rdname NGivenS
qNGivenS <- function(p, S, pars, mod) {
    switch(mod, 'fish' = qnorm(p, .fishMean(pars) * S, sqrt(.fishVar(pars) * S)),
           'plnorm' = 0)
}

# ----
# functions for the means and variances of different SAD models

.fishMean <- function(beta) {
    p <- exp(-beta)

    -1 / log(1 - p) * p / (1 - p)
}

.fishVar <- function(beta) {
    p <- exp(-beta)

    -(p^2 + p * log(1 - p)) / ((1 - p)^2 * (log(1 - p))^2)
}


