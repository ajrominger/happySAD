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
#' @references Levin, B. (1981). A representation for multinomial cumulative distribution functions.
#' The Annals of Statistics, 1123-1126.
#'
#' @export

pmultinom <- function(q, size, prob , s = size) {
    if(size > 0) {
        # s <- ifelse(size < 12, size, 12)
        # s <- size
        prob <- prob / sum(prob)
        # browser()

        o <- prod(ppois(q, s * prob)) / dpois(size, s) *
            dsumuptpois(size, s * prob, q)
        # o[o > 1] <- 1
        # o[o < 0] <- .Machine$double.eps

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
# @references Levin, B. (1981). A representation for multinomial cumulative distribution functions.
# The Annals of Statistics, 1123-1126.
#'
#' @export

# dsumuptpois <- function(N, lambda, m) {
#     fun <- DiscreteDistribution(0:m[1], duptpois(0:m[1], lambda[1], m[1]))
#
#     for(i in 2:length(lambda)) {
#         fun <- fun + DiscreteDistribution(0:m[i], duptpois(0:m[i], lambda[i], m[i]))
#     }
#
#     return(fun@d(N))
# }

dsumuptpois <- function(N, lambda, m) {
    ff <- vector('list', length(lambda))
    ff[[1]] <- duptpois(0:ifelse(m[1] < N, m[1], N), lambda[1], m[1])

    for(i in 2:length(lambda)) {
        ff[[i]] <- cladoRcpp::rcpp_convolve(ff[[i - 1]],
                                            duptpois(0:ifelse(m[i] < N, m[i], N), lambda[i], m[i]))
    }
    # browser()

    o <- ff[[length(lambda)]][N + 1]
    o[is.na(o)] <- 0

    return(o)
}

# dsumuptpois <- function(N, lambda, m) {
#     ff <- vector('list', length(lambda))
#     ff[[1]] <- duptpois(0:m[1], lambda[1], m[1])
#
#     for(i in 2:length(lambda)) {
#         ff[[i]] <- convolve(ff[[i - 1]], rev(duptpois(0:m[i], lambda[i], m[i])), type = 'o')
#     }
#     # browser()
#
#     o <- ff[[length(lambda)]][N + 1]
#     o[is.na(o)] <- 0
#
#     return(o)
# }

# dsumuptpois <- function(N, lambda, m) {
#     # number of random variables
#     t <- length(m)
#
#     # vector of means and variances for each r.v.
#     mu <- lambda * (1 - dpois(m, lambda) / ppois(m, lambda))
#     sig2 <- mu - (m - mu) * (lambda - mu)
#
#     # function for Bruce Levin's $m^r$
#     msupR <- function(r) {
#         sapply(m, function(mi) prod(mi - 0:(r + 1)))
#     }
#
#     # vector of Bruce Levin's $\mu_{(r)}$
#     mus <- matrix(0, nrow = 4, ncol = t)
#     mus[1, ] <- mu
#     for(r in 2:4) {
#         mus[r, ] <- lambda * mus[r - 1, ] - msupR(r - 1) * (lambda - mu)
#     }
#
#     # third and forth moments
#     mu3 <- mus[3, ] + mus[2, ] * (3 - 3*mu) + (mu - 3*mu^2 + 2*mu^3)
#     mu4 <- mus[4, ] + mus[3, ] * (6 - 4*mu) + mus[2, ] * (7 - 12*mu + 6*mu^2) +
#         (mu - 4*mu^2 + 6*mu^3 - 3*mu^4)
#
#     # coeffs of skewness and excess
#     gamma1 <- 1/sqrt(t) * (sum(mu3) / t) / ((sum(sig2) / t)^(3/2))
#     gamma2 <- 1/t * (sum(mu4 - 3 * sig2^2) / t) / ((sum(sig2) / t)^2)
#
#     # Bruce Levin's $f(x)$
#     f <- function(x) {
#         (exp(-(x^2) / 2) / sqrt(2*pi)) * (1 + gamma1/6 * (x^3 - 3*x) +
#                                               gamma2/24 * (x^4 - 6*x^2 + 3) +
#                                               (gamma1^2)/72 * (x^6 - 15*x^4 + 45*x^2 - 15))
#     }
#
#     return(f((N - sum(mu)) / sqrt(sum(sig2))) / sqrt(sum(sig2)))
# }

#' @title Upper truncated Poisson random variables
#'
# @description
#'
# @details
#'
#' @param x vector of integer values at which to calculate probabilities
#' @param lambda mean of the upper truncted Poisson
#' @param m upper (inclusive) limit
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @references Levin, B. (1981). A representation for multinomial cumulative distribution functions.
#' The Annals of Statistics, 1123-1126.
#'
#' @export

duptpois <- function(x, lambda, m) {
    o <- dpois(x, lambda) / ppois(m, lambda)
    o[x > m] <- 0

    return(o)
}


# conv <- function(dfun1, dfun2, m) {
#     M <- sum(m)
#     mm <- 0:M
#
#     function(x) {
#         return(colSums(outer(mm, x, function(m, k) dfun1(m) * dfun2(k - m))))
#     }
# }



# mydsumuptpois <- function(N, lambda, m) {
#     funlist <- vector('list', length(lambda))
#
#     funlist[[1]] <- function(x) duptpois(x, lambda[1], m[1])
#
#     for(i in 2:length(lambda)) {
#         miless1 <- sum(m[1:(i - 1)])
#         funlist[[i]] <- conv(funlist[[i - 1]], function(x) duptpois(x, lambda[i], m[i]), m = c(miless1, m[i]))
#     }
#
#     return(funlist[[length(lambda)]](N))
# }
#





# foo <- conv(conv(conv(function(x) dpois(x, 1),
#                       function(x) dpois(x, 3), 100),
#                  function(x) dpois(x, 5), 100),
#             function(x) dpois(x, 7), 100)
# foo(0:10) - dpois(0:10, 1 + 3 + 5 + 7)
#
# distrdsumuptpois(0:10, c(1, 2, 5, 7), rep(100, 4))
#

#
#
# bla <- conv(conv(function(x) dpois(x, 1), function(x) dpois(x, 2), 1000),
#             function(x) dpois(x, 3), 1000)
# plot(0:10, sapply(0:10, bla))
# lines(0:10, dpois(0:10, 6))

