x <- rfish(1000, 0.01)
xsad <- sad(x, model = 'fish', keepData = TRUE)
# plot(xsad, ptype = 'rad', log = 'y')

N <- c(0, round(exp(seq(log(1), log(0.5*max(x)), length.out = 10))))

sampThr <- lapply(N, function(n) {
    foo <- dhyper(n, 1:(10*max(x)), sum(x) - (1:(10*max(x))), round(0.5*sum(x)))
    thr <- sum((foo * dfish(1:(10*max(x)), xsad$MLE))[is.finite(foo)])
    
    obs <- mean(replicate(100, {
        temp <- rep(1:length(x), x)
        s <- rowSums(outer(1:length(x), sample(temp, round(0.5*sum(x))), '=='))
        mean(s == n)
    }))
    
    return(c(thr = thr, obs = obs))
})
sampThr <- do.call(rbind, sampThr)

N <- N[sampThr[, 2] > 0]
sampThr <- sampThr[sampThr[, 2] > 0, ]

plot(N, sampThr[, 1], ylim = range(sampThr),
     type = 'b', pch = 16, cex = 0.5, xlim = c(0, N[5]))
points(N, sampThr[, 2], type = 'b', cex = 0.6, pch = 16, col = 'red')

plot(sampThr)
abline(0, 1)
