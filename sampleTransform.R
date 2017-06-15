x <- rfish(1000, 0.001)
xsad <- sad(x, model = 'fish', keepData = TRUE)
plot(xsad, ptype = 'rad', log = 'y')

N <- round(exp(seq(log(1), log(0.5*max(x)), length.out = 50)))

sampThr <- lapply(N, function(n) {
    foo <- dhyper(n, 1:(10*max(x)), sum(x) - (1:(10*max(x))), round(0.5*sum(x)))
    thr <- sum((foo * dfish(1:(10*max(x)), xsad$MLE))[is.finite(foo)])
    
    obs <- mean(replicate(100, {
        s <- sample.sad(xsad, round(0.5*sum(x)))
        mean(s == n)
    }))
    
    return(c(thr = thr, obs = obs))
})
sampThr <- do.call(rbind, sampThr)
N <- N[sampThr[, 2] > 0]
sampThr <- sampThr[sampThr[, 2] > 0, ]

plot(N, sampThr[, 1], ylim = range(sampThr), log = 'x',
     type = 'b', pch = 16, cex = 0.5)
points(N, sampThr[, 2], type = 'b', cex = 0.6, pch = 16, col = 'red')

plot(sampThr)
abline(0, 1)
