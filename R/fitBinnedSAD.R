fitBinnedSAD <- function(x, model, keepData = TRUE) {
    init <- switch(model,
                   'fish' = 0.01,
                   'tnegb' = c(1, 1),
                   'plnorm' = c(4, 10),
                   'stick' = 0.8)
    dfun <- switch(model,
                   'fish' = dfish,
                   'tnegb' = dtnegb,
                   'plnorm' = dplnorm,
                   'stick' = dstick)

    dat <- .aggX(x)

    fun2opt <- function(pars) {
        pp <- do.call(dfun, c(list(dat$x), as.list(pars)))
        multip <- as.numeric(tapply(pp, dat$b, sum))

        # cat(pars, range(multip), '', sep = '***')
        if(all(multip < .Machine$double.eps)) {
            o <- 1e+06
        } else {
            o <- - dmultinom(dat$binx, sum(dat$binx), multip, log = TRUE)
        }

        o[o > 1e+06] <- 1e+06
        # cat(o, '\n')

        return(o)
    }

    # browser()

    if(model == 'plnorm') {
        lor <- c(-10, 0.0001)
    } else if(model == 'tnegb') {
        lor <- c(0.0001, 0.0001)
    } else {
        lor <- 0.000001
    }

    if(model == 'stick') {
        upr <- 1
    } else if(model == 'fish') {
        upr <- 2
    } else {
        upr <- c(100, 100)
    }

    # browser()
    optim(init, fun2opt, method = ifelse(length(init) == 1, 'Brent', 'L-BFGS-B'),
          lower = lor, upper = upr)
    # optim(init, fun2opt)
}

.aggX <- function(x) {
    binx <- floor(log(x, 2))
    xx <- 1:(2^(max(binx) + 1) - 1)
    bb <- floor(log(xx, 2))

    binx <- factor(binx, levels = 0:max(binx))

    return(list(binx = as.integer(table(binx)), x = xx, b = bb))
}

# .aggProbs <- function(x, b, dfun) {
#     pp <- dfun(x)
#
#     return(tapply(pp, b, sum))
# }
#
# x <- rtnegb(200, 5, 0.1)
# fitBinnedSAD(x, 'tnegb')[c(1, 2, 4)]
#
#
#
# fitSAD(moths, 'fish')
#
# x <- rfish(200, 0.01)
# fitSAD(x, 'fish')
# fitBinnedSAD(x, 'fish')[c(1, 2, 4)]
# fitBinFish(moths, 0.1)
# lines(fitGambin(x))
# fit_abundances(x)
