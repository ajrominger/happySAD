## ======================================================================
## simulation that generates random SAD data under different conditions:
##   1) different numbers of species
##   2) different subsampling sizes
##   3) different models
##   4) different parameter values
## and then tests:
##   1) AIC methods
##   2) likelihood goodness of fit methods
##   3) mse goodness of fit methods on both rank and cumulative functions
## ======================================================================

setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')


## =================
## simulation set up
## =================

## number of species
nspp <- 10^seq(1.75, 2.75, length=5) - 10^seq(1.75, 2.75, length=5) %% 50

## proportion sampled
prop <- 10^seq(-0.6, 0, length=5) - 10^seq(-0.6, 0, length=5) %% 0.05

## SAD parameters
sad.par <- list(fish=10^seq(-1, -3, length=4),
                plnorm=list(c(0, 1), c(3, 1), c(1, 2), c(0, 2.5)),
                stick=seq(0.5, 1/50, length=4),
                tnegb=list(c(1, 1), c(8, 1), c(20, 1.5), c(100, 1.5)),
                tpois=seq(5, 50, length=4))

sad.rfun <- lapply(names(sad.par), function(f) {
    lapply(1:4, function(p) {
        if(f %in% c('plnorm', 'tnegb')) {
            return(function(n) get(sprintf('r%s', f))(n, sad.par[[f]][[p]][1], sad.par[[f]][[p]][2]))
        } else {
            return(function(n) get(sprintf('r%s', f))(n, sad.par[[f]][p]))
        }
    })
})

names(sad.rfun) <- names(sad.par)

rapply(sad.rfun, function(f) f(3), how='unlist')

## ==============================
## function to run one simulation
## ==============================

simSAD <- function(rfuns, nspp, prop, nrep=1000) {
    rapply(rfuns, function(f) {
        ## simulate data
        dat <- sample.sad(f(nspp), prob=prop)
        
        ## fit SAD models
        fit <- fitSAD(dat, keepData=FALSE)
        
        ## apply over fitted SAD objects
        fit.stats <- sapply(fit, function(x) {
            ## extract needed functions and n from fitted SAD
            n <- x$nobs
            newrfun <- getrfun(x)
            newdfun <- getdfun(x)
            predRank <- sad2Rank(x)
            predCDF <- getpfun(x)(unique(dat))
            
            ## helper function to extract needed stats from simulation
            getStats <- function(r) {
                r <- sort(r, decreasing=TRUE)
                obsCDF <- .ecdf(r)[, 2]
                
                radE <- r - predRank
                radElog <- log(r) - log(predRank)
                cdfE <- obsCDF - predCDF
                cdfElog <- log(obsCDF) - log(predCDF)
                
                c(ll=sum(newdfun(r, log=TRUE)), 
                  radMSE=mean(radE^2), radMSElog=mean(radElog^2), 
                  radMSErel=mean((radE/r)^2), radMSElogrel=mean((radElog/log(r))^2),
                  cdfMSE=mean(cdfE^2), cdfMSElog=mean(cdfElog^2), 
                  cdfMSErel=mean((cdfE/r)^2), cdfMSElogrel=mean((cdfElog/log(obsCDF))^2))
            }
            
            ## simulate for z-scores
            z <- replicate(nrep, getStats(newrfun(x$nobs)))
            z <- rbind(getStats(dat), z)
            
            zOut <- apply(z, 2, function(s) {
                Z <- ((s[1] - mean(s))/sd(s))^2
                P <- sum(Z > ((s - mean(s))/sd(s))^2)/nrep
                return(c(z=Z, p=P))
            })
            
            ## return AIC and summarized z-scores
            return(c(aic=AIC(x), as.vector(zOut)))
        })
        
        return(fit.stats)
    }) 
}



## ======================
## save simulation output
## ======================
save(nspp, prop, sad.par, file='sim_out.RData')
