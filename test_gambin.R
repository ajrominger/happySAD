setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(gambin)
library(parallel)

## funciton to calculate probabilities for binned fisher log series from raw data
dBinFish <- function(x, beta, log=FALSE) {
    maxBinX <- (2^(1+floor(log(max(x), 2)))-1)
    probs <- dfish(1:maxBinX, beta)
    binProbs <- tapply(probs, floor(log(1:maxBinX, 2)), sum)
    binProbs <- binProbs/sum(binProbs)
    newx <- create_octaves(x)[, 2]
    
    if(log) {
        return(log(binProbs) * newx)
    } else {
        return(exp(newx * log(binProbs)))
    }
}

## function to 'unbin' data
unbin <- function(bins, freq) {
    rep(2^bins, round(freq))
}

## fit the binned fisher log series
fitBinFish <- function(x, init) {
    fun <- function(b) -sum(dBinFish(x, b, log=TRUE))
    out <- optimize(fun, interval=c(10^-5, 3*init))
    out <- c(out[[1]], out[[2]])
    names(out) <- c('MLE', 'll')
    out[2] <- -out[2]
    
    return(out)
}


## set-up simulation
bs <- 10^seq(-3.5, 0, length=5)
ns <- 10^seq(1.95, 2.7, length=4)
ns <- ns - ns %% 20
nsim <- 6

gb.sim <- lapply(ns, function(n) {
    b.out <- lapply(bs, function(b) {
        sim.out <- mclapply(1:nsim, mc.cores=6, FUN=function(i) {
            ## simulate from fisher
            r <- rfish(n, b)
            r <- r[is.finite(r)]
            
            ## fit gambin
            gbfit <- fitGambin(r)
            
            gbout <- c(gbfit$Alpha, gbfit$logLik)
            names(gbout) <- c('MLE', 'll')
            
            ## fit fisher to raw
            ffit <- sad(r, 'fish')
            fout <- c(ffit$MLE, ffit$ll)
            names(fout) <- c('MLE', 'll')
            bfrout <- fout
            bfrout[2] <- sum(dBinFish(r, fout[1], log=TRUE))
            
            ## fit tnegb to raw
            tnbfit <- sad(r, 'tnegb')
            tnbout <- c(tnbfit$MLE[2], tnbfit$ll)
            names(tnbout) <- c('MLE', 'll')
            
            ## fit binned fisher
            bfout <- fitBinFish(r, init=ffit$MLE)
            
            out <- cbind(iter=i, nspp=n, truPar=b, rbind(gambin=gbout, fish=fout, rawBinFish=bfrout, binFish=bfout, tnegb=tnbout))
            return(out)
        })
        
        do.call(rbind, sim.out)
    })
    
    do.call(rbind, b.out)
})

gb.sim <- do.call(rbind, gb.sim)

save(gb.sim, file='gambin_sim.RData')