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
#' @return a vector of probabilities equal in length to \code{q}
#'
#' @references Levin, B. (1981). A representation for multinomial cumulative distribution functions.
#' The Annals of Statistics, 1123-1126.
#'
#' @export

pmultinom <- function(q, size, prob, s = size) {
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
#' @param N integer (potentially a vector), the value of summed upper truncated Poisson r.v.'s
#' @param lambda vector of means for each upper truncted Poisson
#' @param m integer vector upper (inclusive) limits
#'
#' @return a vector of probabilities equal in length to \code{N}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

dsumuptpois <- function(N, lambda, m) {
    ff <- vector('list', length(lambda))
    ff[[1]] <- duptpois(0:ifelse(m[1] < N, m[1], N), lambda[1], m[1])

    for(i in 2:length(lambda)) {
        ff[[i]] <- cladoRcpp::rcpp_convolve(ff[[i - 1]],
                                            duptpois(0:ifelse(m[i] < N, m[i], N),
                                                     lambda[i], m[i]))
    }

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
#'
#' @return a vector of probabilities equal in length to \code{x}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

duptpois <- function(x, lambda, m) {
    o <- dpois(x, lambda) / ppois(m, lambda)
    o[x > m] <- 0

    return(o)
}
