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
                                            duptpois(0:ifelse(m[i] < N, m[i], N),
                                                     lambda[i], m[i]))
    }
    # browser()

    o <- ff[[length(lambda)]][N + 1]
    o[is.na(o)] <- 0

    return(o)
}


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

