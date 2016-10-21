library(pika)
library(socorro)

setwd('~/Dropbox/Research/happySAD')


## ============================================================
## functions to calculate different summary statistics for SADs
## ============================================================

chisq.sad <- function(x, type = c('cdf', 'rad')) {
    type <- match.arg(type, c('cdf', 'rad'))
    
    if(type == 'cdf') {
        obs <- simpECDF(x$data)
        O <- obs[, 2]
        E <- getpfun(x)(obs[, 1])
    } else {
        O <- sad2Rank(x)
        E <- sort(x$data, TRUE)
    }
    
    return(sum((O - E)^2 / E))
}

ks.sad <- function(x) {
    obs <- simpECDF(x$data)
    thr <- getpfun(x)(obs[, 1])
    return(max(abs(obs[, 2] - thr)))
}

## likelihood ratio assuming null model is geometic (='broken stick')
lr.sad <- function(x) {
    mod0 <- sad(x$data, 'stick')
    return(x$ll - mod0$ll)
}


## =================================================================
## simulate SADs from different generating models, fitting different
## models to each that have variable parametric agreement
## =================================================================

## generating and fitted models
fg <- list(sad(model = 'tnegb', par = c(2, 0.01)), 
           sad(model = 'tnegb', par = c(6, 1)))
ff <- list(sad(model = 'tnegb', par = c(2, 0.01)), 
           sad(model = 'tnegb', par = c(6, 1)),
           sad(model = 'tnegb', par = c(10, 3)))

## different numbers of species to simulate
S <- round(exp(seq(log(10), log(200), length = 10)) / 5) * 5

## number of data sets to simulate
nrep <- 100

## number of iterations for z-val calculations
nz <- 500

out <- lapply(fg, function(g) {
    ## all data possibly needed (one column per replicate)
    X <- matrix(getrfun(g)(max(S) * nrep), ncol = nrep)
    
    ## loop over fitted models
    outF <- lapply(ff, function(f) {
        ## all possibly needed data simulated from the fitted model
        fdat <- matrix(getrfun(f)(max(S) * nz), ncol = nz)
        
        ## loop over species numbers
        outS <- lapply(S, function(s) {
            ## calculate statistics for fitted model given data from generating model
            singleStat <- sapply(1:nrep, function(i) {
                fm <- f
                fm$data <- X[1:s, i]
                fm$ll <- sum(getdfun(fm)(fm$data, log = TRUE))
                fm$nobs <- s
                
                return(c(ks = ks.sad(fm),
                         ll = logLik(fm),
                         chi2.cdf = chisq.sad(fm, type = 'cdf'),
                         chi2.rad = chisq.sad(fm, type = 'rad'),
                         mse.cdf = mse(fm, type = 'c'),
                         mse.rad = mse(fm, type = 'r'),
                         lr = lr.sad(fm)))
            })
            
            ## calculate statistics for fitted model given data simulated from fitted model
            ## (i.e. for z-value calculations)
            zSim <- sapply(1:nz, function(i) {
                fm <- f
                fm$data <- fdat[1:s, i]
                fm$ll <- sum(getdfun(fm)(fm$data, log = TRUE))
                fm$nobs <- s
                
                return(c(ks = ks.sad(fm),
                         ll = logLik(fm),
                         chi2.cdf = chisq.sad(fm, type = 'cdf'),
                         chi2.rad = chisq.sad(fm, type = 'rad'),
                         mse.cdf = mse(fm, type = 'c'),
                         mse.rad = mse(fm, type = 'r'),
                         lr = lr.sad(fm)))
            })
            
            ## calculate z-val
            zMean <- apply(zSim, 1, mean)
            zSD <- apply(zSim, 1, sd)
            zVal <- ((singleStat - zMean) / zSD)^2
            
            vals <- t(rbind(singleStat, zVal))
            colnames(vals) <- paste(c(rep('', ncol(vals)/2), rep('z_', ncol(vals)/2)), 
                                   colnames(vals), sep = '')
            
            return(vals)
        })
        
        outS <- data.frame(ff = paste(f$MLE, collapse = ','), S = rep(S, each = nrep), do.call(rbind, outS))
        return(outS)
        
    })
    
    outF <- cbind(fg = paste(g$MLE, collapse = ','), do.call(rbind, outF))
    return(outF)
})

out <- do.call(rbind, out)

write.csv(out, file = 'sim_goodnessFitness.csv', row.names = FALSE)
