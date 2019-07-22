dNGivenS <- function(N, S, pars, mod) {
    switch(mod, 'fish' = dnorm(N, .fishMean(pars) * S, sqrt(.fishVar(pars) * S)),
           'plnorm' = 0)
}

pNGivenS <- function(N, S, pars, mod) {
    switch(mod, 'fish' = pnorm(N, .fishMean(pars) * S, sqrt(.fishVar(pars) * S)),
           'plnorm' = 0)
}

qNGivenS <- function(p, S, pars, mod) {
    switch(mod, 'fish' = qnorm(p, .fishMean(pars) * S, sqrt(.fishVar(pars) * S)),
           'plnorm' = 0)
}

.fishMean <- function(beta) {
    p <- exp(-beta)

    -1 / log(1 - p) * p / (1 - p)
}

.fishVar <- function(beta) {
    p <- exp(-beta)

    -(p^2 + p * log(1 - p)) / ((1 - p)^2 * (log(1 - p))^2)
}


dNGivenS(round(seq(100, 2000, length.out = 10)), 20, 0.1, 'fish')
qNGivenS(0.5, 100, 0.1, 'fish')


dSGivenN <- function(S, N, pars, mod) {
    maxS <- round(N / (0.1 * .fishMean(pars)))
    Z <- sum(dNGivenS(N, 1:maxS, pars, mod))
    o <- dNGivenS(N, S, pars, mod) / Z

    o[S >= maxS] <- 0

    return(o)
}
