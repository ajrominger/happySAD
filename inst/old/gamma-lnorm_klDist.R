library(plyr)
library(parallel)
library(viridis)

dDklGamma <- function(x, shape, scale, meanlog, sdlog) {
    dgamma(x, shape, scale = scale) * (dgamma(x, shape, scale = scale, log = TRUE) - 
                                           dlnorm(x, meanlog, sdlog, log = TRUE))
}

DklGamma <- function(shape, scale, meanlog, sdlog) {
    integrate(function(x) dDklGamma(x, shape, scale = scale, meanlog, sdlog), 
              lower = 0, upper = Inf, 
              rel.tol = .Machine$double.eps^0.25)
}


parMin <- 0.1
parMax <- 8
badOptim <- 666666666


parRange <- seq(parMin, parMax, length.out = 5)
par(mfcol = c(length(parRange), length(parRange)), mar = rep(1, 4))
plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '')
for(i in 1:length(parRange)) {
    for(j in 1:length(parRange)) {
        par(mfg = c(length(parRange) - j + 1, i))
        curve(dgamma(x, shape = parRange[i], scale = parRange[j]), from = 0.1, to = 100, 
              ylim = c(0.00001, 0.2), log = 'x')
    }
}



parRange <- seq(parMin, parMax, length.out = 20)
pars <- expand.grid(shape = parRange, scale = parRange)
gammaAsLnorm <- mclapply(1:nrow(pars), mc.cores = 6, function(i) {
    out <- optim(par = c(1, 1), method = 'Nelder-Mead',
                 fn = function(lnormPars) {
                     x <- try(DklGamma(pars[i, 1], pars[i, 2], 
                                       meanlog = lnormPars[1], sdlog = lnormPars[2])$value, 
                              silent = TRUE)
                     if('try-error' %in% class(x)) {
                         return(badOptim)
                     } else {
                         return(x)
                     }
                 })
    if(out$value == badOptim) {
        out <- optim(par = c(1, 1), method = 'SANN',
                     fn = function(lnormPars) {
                         x <- try(DklGamma(pars[i, 1], pars[i, 2], 
                                           meanlog = lnormPars[1], sdlog = lnormPars[2])$value, 
                                  silent = TRUE)
                         if('try-error' %in% class(x)) {
                             return(badOptim)
                         } else {
                             return(x)
                         }
                     })
    }
    
    return(unlist(out[1:2]))
})

gammaAsLnorm <- as.data.frame(do.call(rbind, gammaAsLnorm))
gammaAsLnorm$value[gammaAsLnorm$value == badOptim] <- NA

par(mfcol = c(1, 1))
image(parRange, parRange, matrix(gammaAsLnorm$value, nrow = length(parRange)), col = viridis(60))
range(gammaAsLnorm$value, na.rm = TRUE)


mu <- seq(-13, 8, length.out = 20)
sig <- seq(0.1, 24, length.out = 20)
MuSig <- expand.grid(mu, sig)
foo <- lapply(1:nrow(MuSig), function(i) {
    DklGamma(shape = parRange[1], scale = parRange[5], meanlog = MuSig[i, 1], sdlog = MuSig[i, 2])
})

