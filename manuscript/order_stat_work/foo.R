ks.sad <- function(x) {
    dat <- x$data
    
    # cumulative density function
    pfun <- getpfun(x)
    
    # cumulative prob observed and from theory
    n <- length(dat)
    pobs <- (n:1) / n
    pthr <- pfun(sort(dat))
    
    # the statisic is the difference between obs and thr
    out <- pthr - pobs
    
    return(max(out, 1 / n - out, na.rm = TRUE))
}


# number of parameters to use
nparam <- 25

# unique sample sizes
nn <- seq(30, 200, length.out = 6)

# number of simulations to run
nsim <- 100

# the generating models
genMods <- rbind(expand.grid(mod = 'rplnorm', 
                             n = nn,
                             p1 = seq(1, 8, length.out = sqrt(nparam)), 
                             p2 = seq(0.1, 4, length.out = sqrt(nparam))), 
                 expand.grid(mod = 'rtnegb', 
                             n = nn, 
                             p1 = seq(1, 8, length.out = sqrt(nparam)), 
                             p2 = seq(0.01, 4, length.out = sqrt(nparam))), 
                 expand.grid(mod = 'rfish', 
                             n = nn,
                             p1 = seq(0.0001, 0.1, length.out = nparam), 
                             p2 = NA))
genMods$mod <- as.character(genMods$mod)


goodness <- function(x) {
    o <- sapply(x, function(thisSAD) {
        c('logLikZ' = logLikZ(thisSAD), 'mseZ' = mseZ(thisSAD, 4, relative = FALSE, log = FALSE), 
          'mse' = mse(thisSAD, relative = FALSE, log = FALSE), 'ks' = ks.sad(thisSAD), 
          'logLik' = logLik(thisSAD), 'AIC' = AIC(thisSAD))
    })
    o <- data.frame(fitted_mod = colnames(o), t(o))
    rownames(o) <- NULL
    return(o)
}

allGood <- parallel::mclapply(1:nrow(genMods), mc.cores = 3, function(i) {
    # repeat nsim times
    sims <- replicate(1, simplify = FALSE, expr = {
        # one realization
        x <- get(genMods[i, 1])(genMods[i, 2], genMods[i, 3], genMods[i, 4])
        
        # calcualte goodness of fit metrics
        o <- goodness(fitSAD(x, c('fish', 'plnorm', 'tnegb'), keepData = TRUE))
        o <- cbind(gen_mod = substring(genMods[i, 1], 2), o)
        
        return(o)
    })
    
    sims <- do.call(rbind, sims)
    return(sims)
})

allGood <- do.call(rbind, allGood)
for(j in 3:8) allGood[, j] <- as.numeric(allGood[, j])

pairs(allGood[, 3:6], log = 'xy')
