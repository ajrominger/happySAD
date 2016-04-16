setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(gambin)

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
bs <- 10^seq(-3.25, 0, length=5)

gb.sim <- lapply(bs, function(b) {
    ## simulate from fisher
    r <- rfish(100, b)
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
    
    ## fit binned fisher
    bfout <- fitBinFish(r, init=ffit$MLE)
    
    return(rbind(gambin=gbout, fish=fout, rawBinFish=bfrout, binFish=bfout))
})

fitBinFish(r)
sad(r, 'fish')




bla <- gambin_exp(7.1, 11, 100)
ubla <- unbin((1:length(bla)) - 1, bla)
ublaf <- sad(ubla, 'fish')

logLik(fitGambin(ubla))
dgambin

plot(dBinFish(2^((1:length(bla)) - 1), ublaf$MLE))
points(dgambin(7.1, 11))
