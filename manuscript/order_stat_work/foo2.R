dorder.sad <- function(r, x, s) {
    # browser()
    n <- s$nobs
    k <- n - r + 1
    dfun <- getdfun(s)
    pfun <- getpfun(s)
    
    jj <- 0:(n - k)
    
    o <- sapply(x, function(xi) {
        sum(choose(n, jj) * (
            (1 - pfun(xi))^jj * 
                pfun(xi)^(n - jj) - 
                (1 - pfun(xi) + dfun(xi))^jj * 
                (pfun(xi) - dfun(xi))^(n - jj)
        ))
    })
    
    o[o < .Machine$double.eps^0.75] <- 0
    
    return(o)
}


dmse <- function(s) {
    # browser()
    rad <- sad2Rank(s)
    
    xx <- 1:(10*max(rad))
    rr <- 1:s$nobs
    
    ords <- sapply(rr, dorder.sad, x = xx, s = s)
    errs <- outer(xx, rad, function(x, y) (x - y)^2)
    
    mean(colSums(ords * errs))
}

foo <- sad(rfish(100, 0.1), 'fish', keepData = TRUE)
z <- mseZ(foo, 999, return.sim = TRUE)
plot(density(log(z$sim)))
log(dmse(foo))




r <- 30
layout(matrix(1:2, nrow = 1))
par(mar = c(2, 2, 0, 0))
bla <- dorder.sad(r, 1:max(foo$data), foo)
plot(bla, 1:max(foo$data), log = 'y', type = 'l', lwd = 2)
plot(foo, ptype = 'rad', log = 'y')
abline(v = r, h = which.max(bla))
