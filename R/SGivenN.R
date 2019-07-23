#' @title Probability of S given N
#'
#' @description Density function for the number of species given the total number of sampled individuals
#'
# @details
#'
#' @param S number of species
#' @param N total abundance
#' @param pars parameter(s) for the SAD model
#' @param mod a string specifying the SAD model
#'
#' @return A numeric vector of length equal to the input
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export
#' @rdname SGivenN

dSGivenN <- function(S, N, pars, mod) {
    maxS <- round(N / (0.1 * .fishMean(pars)))
    Z <- sum(dNGivenS(N, 1:maxS, pars, mod))
    o <- dNGivenS(N, S, pars, mod) / Z

    o[S >= maxS] <- 0

    return(o)
}
