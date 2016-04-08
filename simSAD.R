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
# devtools::load_all('../pika')
# install.packages('~/Dropbox/Research/pika', repos=NULL, type='source')
devtools::install_github('ajrominger/pika')
# library(pika)
library(parallel)


## =================
## simulation set up
## =================

## number of species
nspp <- 10^seq(1.95, 2.7, length=4)
nspp <- nspp - nspp %% 20

## proportion sampled
prop <- 10^seq(-0.4, 0, length=4)
prop <- prop - prop %% 0.05

## SAD parameters
sad.par <- list(fish=10^seq(-1.25, -2, length=4),
                plnorm=list(c(2, 0.5), c(3, 0.5), c(0, 1.5), c(1, 1.5)),
                stick=10^seq(-0.7, -1.7, length=4),
                tnegb=list(c(20, 4), c(8, 1.5), c(8, 0.5), c(100, 3))
                # tpois=seq(5, 50, length=4)
                )

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




## ==============================
## function to run one simulation
## ==============================

simSAD <- function(rfuns, nspp, prop, nrep=1000) {
    rapply(rfuns, how='replace', f=function(f) {
        ## simulate data
        dat <- sample.sad(f(nspp), prob=prop)
        
        ## fit SAD models
        fit <- fitSAD(dat, keepData=FALSE)
        
        ## names for output
        outNames <- c('aic', 'z_ll', 'p_ll', 'z_radMSE', 'p_radMSE', 'z_radMSElog', 'p_radMSElog', 
                      'z_radMSErel', 'p_radMSErel', 'z_cdfMSE', 'p_cdfMSE', 'z_cdfMSElog', 
                      'p_cdfMSElog', 'z_cdfMSErel', 'p_cdfMSErel')
        
        ## apply over fitted SAD objects
        fit.stats <- sapply(fit, function(x) {
            ## only if model is tpois or stick will we do z-score stuff
            if(x$model %in% c('tpois', 'stick')) {
                ## extract needed functions and n from fitted SAD
                n <- x$nobs
                newrfun <- getrfun(x)
                newdfun <- getdfun(x)
                predRank <- sad2Rank(x, x$nobs)
                predCDF <- getpfun(x)
                
                ## helper function to extract needed stats from simulation
                getStats <- function(r) {
                    r <- sort(r, decreasing=TRUE)
                    if(!is.finite(mean(r)) | !is.finite(var(r))) browser()
                    obsCDF <- .ecdf(r)
                    
                    radE <- r - predRank
                    radElog <- log(r) - log(predRank)
                    cdfE <- obsCDF[, 2] - predCDF(obsCDF[, 1])
                    cdfElog <- log(obsCDF[, 2]) - log(predCDF(obsCDF[, 1]))
                    
                    c(ll=sum(newdfun(r, log=TRUE)), 
                      radMSE=mean(radE^2), radMSElog=mean(radElog^2), 
                      radMSErel=mean((radE/r)^2), # radMSElogrel=mean((radElog/log(r))^2),
                      cdfMSE=mean(cdfE^2), cdfMSElog=mean(cdfElog^2), 
                      cdfMSErel=mean((cdfE/predCDF(obsCDF[, 1]))^2) #, cdfMSElogrel=mean((cdfElog/log(predCDF(obsCDF[, 1])))^2)
                    )
                }
                
                ## simulate for z-scores
                z <- replicate(nrep, getStats(newrfun(x$nobs)))
                z <- cbind(getStats(dat), z)
                
                zOut <- apply(z, 1, function(s) {
                    Z <- ((s[1] - mean(s))/sd(s))^2
                    P <- sum(Z <= ((s - mean(s))/sd(s))^2)/(nrep+1)
                    return(c(z=Z, p=P))
                })
                
                ## return AIC and summarized z-scores
                # outNames <- c('aic', as.vector(outer(rownames(zOut), colnames(zOut), paste, sep='_')))
                out <- c(AIC(x), as.vector(zOut))
                
            } else {
                out <- c(AIC(x), rep(NA, length(outNames)-1))
            }
            names(out) <- outNames
            return(out)
        })
        
        return(fit.stats)
    }) 
}

## ==============
## run simulation
## ==============

nsim <- 500
sim.out <- mclapply(1:nsim, mc.cores=3, FUN=function(n) {
    print(n)
    lapply(nspp, function(ns) {
        lapply(prop, function(p) {
            simSAD(sad.rfun, ns, p, 1000)
        })
    })
})

## save simulation output
save(sim.out, nspp, prop, sad.par, sad.rfun, file='sim_out.RData')
