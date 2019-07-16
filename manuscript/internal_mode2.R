library(pika)
bb <- 0.01
nprop <- seq(1, 1, length.out = 1)
nsp <- 100

foo <- lapply(bb, function(b) {
    x <- rfish(1000000, b)
    # print(any(!is.finite(x)))
    x[!is.finite(x)] <- max(x[is.finite(x)])
    sadList <- split(x, ceiling(seq_along(x) / nsp))
    
    samps <- sapply(nprop, function(np) {
        intMode <- sapply(sadList, function(thisOne) {
            # n <- round(np * sum(thisOne))
            # thisOne <- sample.sad(sad(thisOne, keepData = TRUE), n, replace = FALSE)
            octCount <- table(floor(log(thisOne, 2)))
            # browser()
            return(c(any(octCount[1] < octCount[-1]), any(octCount[1] < octCount[2:6])))
        })
        # browser()
        return(rowMeans(intMode))
    })
    
    return(samps)
})

foo

dintModeGivenS <- function(S, b) {
    # browser()
    nn <- 1:(2^13 - 1)
    octs <- floor(log(nn, 2))
    probs <- dfish(nn, b)
    octProbs <- tapply(probs, octs, sum)
    octProbs <- c(octProbs, '13' = 0)
    octProbs <- octProbs / sum(octProbs)
    
    y <- 1:S
    
    o <- sapply(2:length(octProbs), function(i) {
        pbinom(floor((y - 1) / 2), y, octProbs[1] / sum(octProbs[c(1, i)])) *
            dbinom(y, S, sum(octProbs[c(1, i)]))
    })
    
    # o1 <- pbinom(floor((y - 1) / 2), y, octProbs[1] / sum(octProbs[c(1, 2)])) *
    #     dbinom(y, S, sum(octProbs[c(1, 2)]))
    # o2 <- pbinom(floor((y - 1) / 2), y, octProbs[1] / sum(octProbs[c(1, 3)])) *
    #     dbinom(y, S, sum(octProbs[c(1, 3)]))
    # 
    # o3 <- pbinom(floor((y - 1) / 3), y, octProbs[1] / sum(octProbs[c(1, 2, 3)])) *
    #     dbinom(y, S, sum(octProbs[c(1, 2, 3)]))
    
    # return(sum(o1) + sum(o2) - sum(o3))
    return(sum(o))
}



foo
dintModeGivenS(nsp, bb)

dintModeGivenS <- function(S, b) {
    nn <- 1:(2^13 - 1)
    octs <- floor(log(nn, 2))
    probs <- rfish(nn, b)
    octProbs <- tapply(probs, octs, sum)
    octProbs <- c(octProbs, '13' = 0)
    octProbs <- octProbs / sum(octProbs)
    p0 <- octProbs[1]
    
    dS0i <- function(S0i, i) {
        dbinom(S0i, S, p0 + octProbs[i])
    }
    
    db0lessi <- function(S0i, i) {
        pbinom(floor(S0i / 2), S0i, p0 / (p0 + octProbs[i]))
    }
    
    y <- 1:S
    o <- sapply(2:length(octProbs), function(i) {
        db0lessi(y, i) * dS0i(y, i)
    })
    
    return(sum(o))
}

dintModeGivenS(nsp, bb)




.dsumtpois <- function(N, lambda, m) {
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





library(distr)
library(socorro)
.duptpoisConv <- function(x, lambda, m) {
    M <- 0:sum(m)
    sapply(x, function(k) {
        sum(.duptpois(M, lambda[1], m[1]) * .duptpois(k - M, lambda[2], m[2]))
    })
}

.duptpois <- function(x, lambda, m, log = FALSE) {
    `%opp%` <- ifelse(log, `-`, `/`)
    
    o <- dpois(x, lambda, log) %opp% ppois(m, lambda, log = log, lower.tail = TRUE)
    o[x > m | x < 0] <- ifelse(log, -Inf, 0)
    
    return(o)
}

la <- c(9, 10)
mm <- c(5, 20)

dfoo1 <- DiscreteDistribution(0:mm[1], .duptpois(0:mm[1], la[1], mm[1]))
dfoo2 <- DiscreteDistribution(0:mm[2], .duptpois(0:mm[2], la[2], mm[2]))
dfooConv <- DiscreteDistribution(0:sum(mm), .duptpoisConv(0:sum(mm), la, mm))

rfoo <- dfoo1@r(10000) + dfoo2@r(10000)
plot(simpECDF(rfoo))
lines(0:25, dfooConv@p(0:25))


pmultinom <- function(q, size, prob, s = size) {
    prob <- prob / sum(prob)
    
    factorial(size) / (s^size * exp(-s)) * prod(ppois(q, s * prob)) * 
        .dsumtpois(size, s * prob, q)
}


pmultinom(rep(3, 12), 12, rep(1, 12), 12)





nn <- 0:200
ll <- seq(1, 5, length.out = 20)
fooThr <- sapply(nn, .dsumtpois, lambda = ll, m = rep(100, length(ll)))
plot(nn, fooThr, ylim = c(0, max(fooThr)))
points(nn, dpois(nn, sum(ll)), col = 'red')



