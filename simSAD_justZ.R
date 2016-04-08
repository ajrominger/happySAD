setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
# install.packages('~/Dropbox/Research/pika', repos=NULL, type='source')
# library(pika)
library(parallel)


## names for output
outNames <- c('aic', 'z_ll', 'p_ll', 'z_radMSE', 'p_radMSE', 'z_radMSElog', 'p_radMSElog', 
              'z_radMSErel', 'p_radMSErel', 'z_cdfMSE', 'p_cdfMSE', 'z_cdfMSElog', 
              'p_cdfMSElog', 'z_cdfMSErel', 'p_cdfMSErel')

nsim <- 500
nrep <- 500

sim.outZ <- mclapply(1:nsim, mc.cores=3, FUN=function(i) {
    print(i)
    
    dat <- rtnegb(100, 8, 5)
    fit <- fitSAD(dat, c('tpois', 'tnegb'))

    
    fit.stats <- sapply(fit, function(x) {
        ## extract needed functions and n from fitted SAD
        n <- x$nobs
        newrfun <- getrfun(x)
        newdfun <- getdfun(x)
        predRank <- sad2Rank(x, x$nobs)
        predCDF <- getpfun(x)
        
        ## helper function to extract needed stats from simulation
        getStats <- function(r) {
            r <- sort(r, decreasing=TRUE)
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
        names(out) <- outNames
        return(out)
    })
})

save(sim.outZ, file='sim_out_Z.RData')