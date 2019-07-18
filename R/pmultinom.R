#' @title Multinomial CDF
#'
#' @description Cumulative probability function of the multinomial distribution
#'
# @details
#'
#' @param q vector of quantiles, k-long
#' @param size total number of objects put into k boxes
#' @param prob vector of probabilities for each of k boxes, internally normalized
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @references Levin, B. (1981). A representation for multinomial cumulative distribution functions. The Annals of Statistics, 1123-1126.
#'
#' @export

pmultinom <- function(q, size, prob #, s = size
                      ) {
    if(size > 0) {
        s <- ifelse(size < 12, size, 12)
        prob <- prob / sum(prob)
        # browser()

        o <- factorial(size) / (s^size * exp(-s)) * prod(ppois(q, s * prob)) *
            dsumuptpois(size, s * prob, q)
        o[o > 1] <- 1
        o[o < 0] <- .Machine$double.eps

        return(o)
    } else {
        return(1)
    }
}

#' @title Sum of upper truncated Poisson random variables
#'
# @description
#'
# @details
#'
#' @param N integer, the value of summed upper truncated Poisson r.v.'s
#' @param lambda vector of means for each upper truncted Poisson
#' @param m integer vector upper (inclusive) limits
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @references Levin, B. (1981). A representation for multinomial cumulative distribution functions. The Annals of Statistics, 1123-1126.
#'
#' @export

dsumuptpois <- function(N, lambda, m) {
    # number of random variables
    t <- length(m)

    # vector of means and variances for each r.v.
    mu <- lambda * (1 - dpois(m, lambda) / ppois(m, lambda))
    sig2 <- mu - (m - mu) * (lambda - mu)

    # function for Bruce Levin's $m^r$
    msupR <- function(r) {
        sapply(m, function(mi) prod(mi - 0:(r + 1)))
    }

    # vector of Bruce Levin's $\mu_{(r)}$
    mus <- matrix(0, nrow = 4, ncol = t)
    mus[1, ] <- mu
    for(r in 2:4) {
        mus[r, ] <- lambda * mus[r - 1, ] - msupR(r - 1) * (lambda - mu)
    }

    # third and forth moments
    mu3 <- mus[3, ] + mus[2, ] * (3 - 3*mu) + (mu - 3*mu^2 + 2*mu^3)
    mu4 <- mus[4, ] + mus[3, ] * (6 - 4*mu) + mus[2, ] * (7 - 12*mu + 6*mu^2) +
        (mu - 4*mu^2 + 6*mu^3 - 3*mu^4)

    # coeffs of skewness and excess
    gamma1 <- 1/sqrt(t) * (sum(mu3) / t) / ((sum(sig2) / t)^(3/2))
    gamma2 <- 1/t * (sum(mu4 - 3 * sig2^2) / t) / ((sum(sig2) / t)^2)

    # Bruce Levin's $f(x)$
    f <- function(x) {
        (exp(-(x^2) / 2) / sqrt(2*pi)) * (1 + gamma1/6 * (x^3 - 3*x) +
                                              gamma2/24 * (x^4 - 6*x^2 + 3) +
                                              (gamma1^2)/72 * (x^6 - 15*x^4 + 45*x^2 - 15))
    }

    return(f((N - sum(mu)) / sqrt(sum(sig2))) / sqrt(sum(sig2)))
}




dsumuptpois <- function(N, lambda, m) {

}

duptpois <- function(x, lambda, m) {
    dpois(x, lambda) / ppois(m, lambda)
}


conv <- function(dfun1, dfun2, m) {
    M <- sum(m)
    mm <- 1:M
    function(x) sum(dfun1(mm) * dfun2(x - mm))
}


