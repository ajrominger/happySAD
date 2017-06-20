## n is the abundance(s) for which a probability is desired
## x is an SAD object
## N is the number of individuals sampled from that SAD
sampTrans <- function(n, x, N, include0 = TRUE, log = FALSE) {
    if(is.null(x$model)) {
        dfun <- function(n) {
            tab <- table(x$data)
            supp <- as.numeric(names(tab))
            prob <- as.numeric(tab) / sum(tab)
            
            prob <- prob[match(n, supp)]
            prob[is.na(prob)] <- 0
            
            return(prob)
        }
    } else {
        dfun <- getdfun(x)
    }
    
    latentSAD <- dfun((1:sum(x$data)))
    
    nloop <- unique(n)
    
    out <- lapply(nloop, function(i) {
        ## dhyper param conversions
        ## x = n; abundance
        ## m = k; total abundance of focal sp in the complete sample
        ## n = J - k; the abundance of all other spp
        ## k = J; total size of sample (sum of all spp abundances)
        trans <- dhyper(i, 1:sum(x$data), sum(x$data) - (1:sum(x$data)), N)
        
        ## rational:
        ## `trans` is the probability that a species with abundance `k` ends up
        ## with abundance `n` in the sample; we integrate over all `k` while 
        ## multiplying the probability of `n` given `k` by the probability of `k`
        ## given by the SAD
        return(sum((latentSAD * trans)[is.finite(trans)]))
    })
    out <- unlist(out)
    out <- out[match(n, nloop)]
    
    if(log) out <- log(out)
    
    if(include0) {
        return(out)
    } else {
        p0 <- sampTrans(0, x, N, include0 = TRUE)
        if(log) {
            out[n == 0] <- -Inf
            return(out - log(1 - p0))
        } else {
            out[n == 0] <- 0
            return(out / (1 - p0))
        }
    }
}

## testing
# plot(1:30, sampTrans(1:30, xsad, round(0.999*sum(xsad$data)), include0 = FALSE), log = 'y')
# points(1:30, getdfun(xsad)(1:30), col = 'red')
# log(sampTrans(0:4, xsad, 1000, include0 = FALSE))
# sampTrans(0:4, xsad, 1000, log = TRUE, include0 = FALSE)
