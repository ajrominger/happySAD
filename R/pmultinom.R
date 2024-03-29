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
#' @references Levin, B. (1981). A representation for multinomial cumulative distribution
#' functions. The Annals of Statistics, 1123-1126.
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

pmultinomapprox <- function(q, size, prob, s = size) {
    if(size > 0) {
        # s <- ifelse(size < 12, size, 12)
        # s <- size
        prob <- prob / sum(prob)
        # browser()

        o <- prod(ppois(q, s * prob)) / dpois(size, s) *
            dsumapprox(size, s * prob, q)
        # o[o > 1] <- 1
        # o[o < 0] <- .Machine$double.eps

        browser()
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
#'
#' @return a vector of probabilities equal in length to \code{N}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

dsumuptpois <- function(N, lambda, m) {
    # browser()
    appx <- dsumapprox(N, lambda, m)
    appx[is.nan(appx) | is.na(appx)] <- 0

    if(appx < 0.01 | appx > 1 - 0.01) {
        appx[appx < 0] <- 0
        appx[appx > 1] <- 1

        return(appx)
    } else {

        ff <- vector('list', length(lambda))
        ff[[1]] <- duptpois(0:ifelse(m[1] < N, m[1], N), lambda[1], m[1])

        # browser()
        for(i in 2:length(lambda)) {
            ff[[i]] <- cladoRcpp::rcpp_convolve(ff[[i - 1]],
                                                duptpois(0:ifelse(m[i] < N, m[i], N),
                                                         lambda[i], m[i]))
        }

        o <- ff[[length(lambda)]][N + 1]

        o[is.na(o)] <- 0

        return(o)
    }
}

dsumapprox <- function(N, lambda, m) {
    mus <- lambda * (1 - dpois(m, lambda) / ppois(m, lambda))
    sig2s <- mus - (m - mus) * (lambda - mus)
    sig2s[sig2s < 0] <- 0

    dnorm(N, sum(mus), sqrt(sum(sig2s)))
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


# S <- 500
# b <- 0.005
# pp <- tapply(dfish(1:(2^13 - 1), b), floor(log(1:(2^13 - 1), 2)), sum)
# poss <- round(S / 13):S
#
#
# exact <- sapply(poss, function(thisS) pmultinom(rep(thisS, 12), S - thisS, pp[-1]))
# appro <- sapply(poss, function(thisS) pmultinomapprox(rep(thisS, 12), S - thisS, pp[-1]))
#
# plot(appro, abs(exact - appro), xlim = c(0, 0.8))
# abline(h = 0, v = 0.01, col = 'red')
