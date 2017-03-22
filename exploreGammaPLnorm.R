library(socorro)

npar <- 5
mu <- seq(0, 4, length = npar)
sig <- seq(0.1, 3.7, length = npar) + 0.4

par(mfrow = rep(npar, 2), mar = rep(0.5, 4), oma = c(3, 3, 0, 0))
for(i in npar:1) {
    for(j in 1:npar) {
        curve(dlnorm(x, mu[j], sig[i]), from = 0.1, to = 500, axes = FALSE, log = 'x')
        box()
        if(i == 1) logAxis(1, expLab = TRUE)
    }
}

s <- seq(0.1, 10, length = npar)
r <- seq(0.1, 10, length = npar)
par(mfrow = rep(npar, 2), mar = rep(0.5, 4), oma = c(3, 3, 0, 0))
for(i in npar:1) {
    for(j in 1:npar) {
        curve(dgamma(x, s[j], rate = r[i]), from = 0.1, to = 500, axes = FALSE, log = 'x')
        box()
        if(i == 1) logAxis(1, expLab = TRUE)
    }
}