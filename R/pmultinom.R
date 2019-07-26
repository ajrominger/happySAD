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

pmultinomApprox <- function(q, size, prob, s = size) {
    if(size > 0) {
        # s <- ifelse(size < 12, size, 12)
        # s <- size
        prob <- prob / sum(prob)
        # browser()

        o <- prod(ppois(q, s * prob)) / dpois(size, s) *
            dsumapprox(size, s * prob, q)
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

dsumuptpoisCut <- function(N, lambda, m) {
    ff <- vector('list', length(lambda))
    nmin <- 0
    nmax <- ifelse(m[1] < N, m[1], N)
    ff[[1]] <- duptpois(nmin:nmax, lambda[1], m[1])

    # browser()
    for(i in 2:length(lambda)) {
        nmax <- ifelse(m[i] < N, m[i], N)
        ff[[i]] <- cladoRcpp::rcpp_convolve(ff[[i - 1]], duptpois(nmin:nmax, lambda[i], m[i]))
        grtr0 <- ff[[i]] > .Machine$double.eps^0.75
        nmin <- nmin + min(which(grtr0)) - 1
        ff[[i]] <- ff[[i]][grtr0]
    }

    o <- ff[[length(lambda)]][N + 1 - nmin]
    browser()

    o[is.na(o)] <- 0

    return(o)
}

dsumuptpois(400, rep(20, 20), rep(28, 20))
dsumuptpoisCut(400, rep(20, 20), rep(28, 20))

dsumapprox <- function(N, lambda, m) {
    mus <- lambda * (1 - dpois(m, lambda) / ppois(m, lambda))
    sig2s <- mus - (m - mus) * (lambda - mus)

    dnorm(N, sum(mus), sqrt(sum(sig2s)))
}

dsumlevin <- function(N, lambda, m) {
    # number of random variables
    t <- length(m)

    # vector of means and variances for each r.v.
    mu <- lambda * (1 - dpois(m, lambda) / ppois(m, lambda)) #CORRECT
    sig2 <- mu - (m - mu) * (lambda - mu) #CORRECT

    # browser()
    # function for Bruce Levin's $m^r$ #CORRECT
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
    gamma1 <- 1/sqrt(t) * mean(mu3) / (mean(sig2)^(3/2)) #CORRECT
    gamma2 <- 1/t * mean(mu4 - 3 * sig2^2) / (mean(sig2)^2) #CORRECT

    # Bruce Levin's $f(x)$ #CORRECT
    f <- function(x) {
        (exp(-(x^2) / 2) / sqrt(2*pi)) * (1 + gamma1/6 * (x^3 - 3*x) +
                                              gamma2/24 * (x^4 - 6*x^2 + 3) +
                                              (gamma1^2)/72 * (x^6 - 15*x^4 + 45*x^2 - 15))
    }

    return(f((N - sum(mu)) / sqrt(sum(sig2))) / sqrt(sum(sig2))) #CORRECT
}

x <- proc.time()
foo <- replicate(100, dsumuptpois(100, seq(5, 10, length.out = 13), rep(20, 13)))
proc.time() - x

x <- proc.time()
foo <- replicate(100, dsumsimp(100, seq(5, 10, length.out = 13), rep(20, 13)))
proc.time() - x



dsumuptpois(50, seq(1, 10, length.out = 12), rep(12, 12))
dsumapprox(50, seq(1, 10, length.out = 12), rep(12, 12))
dsumlevin(50, seq(1, 10, length.out = 12), rep(12, 12))

dsumuptpois(20, rep(2, 12), rep(12, 12))
dsumapprox(20, rep(2, 12), rep(12, 12))
dsumlevin(20, rep(2, 12), rep(12, 12))

dsumuptpois(100, 0.8 * seq(4, 12, length.out = 13), 2 + round(seq(4, 12, length.out = 13)))
dsumapprox(100, 0.8 * seq(4, 12, length.out = 13), 2 + round(seq(4, 12, length.out = 13)))




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
