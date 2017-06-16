sampTrans <- function(n, x, N, include0 = TRUE) {
    dfun <- getdfun(x)
    latentSAD <- dfun((1:sum(x$data)))
    
    ## dhyper param conversions
    ## x = n; abundance
    ## m = k; total abundance of focal sp in the complete sample
    ## n = J - k; the abundance of all other spp
    ## k = J; total size of sample (sum of all spp abundances)
    
    out <- lapply(n, function(i) {
        trans <- dhyper(i, 1:sum(x$data), sum(x$data) - (1:sum(x$data)), N)
        
        ## rational:
        ## `trans` is the probability that a species with abundance `k` ends up
        ## with abundance `n` in the sample; we integrate over all `k` while 
        ## multiplying the probability of `n` given `k` by the probability of `k`
        ## given by the SAD
        return(sum((latentSAD * trans)[is.finite(trans)]))
    })
    out <- unlist(out)
    
    if(include0) {
        return(out)
    } else {
        out[n == 0] <- 0
        p0 <- sampTrans(0, x, N, include0 = TRUE)
        return(out / (1 - p0))
    }
}

## testing
# plot(1:30, sampTrans(1:30, xsad, round(0.999*sum(xsad$data)), include0 = FALSE), log = 'y')
# points(1:30, getdfun(xsad)(1:30), col = 'red')
