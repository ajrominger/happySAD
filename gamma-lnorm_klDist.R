library(plyr)
library(parallel)
library(viridis)

dDklGamma <- function(x, shape, rate, meanlog, sdlog) {
    dgamma(x, shape, rate) * (dgamma(x, shape, rate, log = TRUE) - dlnorm(x, meanlog, sdlog, log = TRUE))
}

DklGamma <- function(shape, rate, meanlog, sdlog) {
    integrate(function(x) dDklGamma(x, shape, rate, meanlog, sdlog), lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.8)
}

parRangeLong <- seq(0.01, 100, length.out = 20)
parRangeShrt <- seq(0.01, 100, length.out = 6)
pars <- expand.grid(shape = parRangeLong, rate = parRangeLong, meanlog = parRangeShrt, sdlog = parRangeShrt)

gammaDkl <- unlist(mclapply(1:nrow(pars), mc.cores = 6, FUN = function(i) {
    out <- try(DklGamma(pars[i, 1], pars[i, 2], pars[i, 3], pars[i, 4])$value, silent = TRUE)
    if('try-error' %in% class(out)) out <- NA
    return(out)
}))

image(matrix(gammaDkl[pars$meanlog == parRangeShrt[4] & pars$sdlog == parRangeShrt[4]], nrow = length(parRangeLong)), 
      col = viridis(60))
contour(matrix(gammaDkl[pars$meanlog == parRangeShrt[4] & pars$sdlog == parRangeShrt[4]], nrow = length(parRangeLong))^0.5)
