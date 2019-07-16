library(parallel)

oneShot <- function(s0, b, pp) {
    x <- sad(rtnegb(s0, 5, 0.01), 'fish', keepData = TRUE)
    N <- sum(x$data)
    n <- round(N * pp)
    
    o <- mclapply(n, mc.cores = 3, mc.preschedule = FALSE, function(ni) {
        newx <- sad(sample.sad(x, ni), 'fish', keepData = TRUE)
        return(c(mle = newx$MLE, ll_z = logLikZ(newx), s0 = newx$nobs))
    })
    o <- do.call(rbind, o)
    
    return(data.frame(n = n, o))
}

itrShot <- function(s0, b, pp) {
    x <- sad(rfish(s0, b), 'fish', keepData = TRUE)
    N <- sum(x$data)
    n <- round(N * pp)
    
    o <- matrix(NA, nrow = length(n), ncol = 4)
    colnames(o) <- c('n', 'mle', 'll_z', 's0')
    
    newx <- x
    o[1, ] <- c(n[1], newx$MLE, logLikZ(newx), newx$nobs)
    
    for(i in 2:length(n)) {
        newx <- sad(sample.sad(newx, n[i]), 'fish', keepData = TRUE)
        o[i, ] <- c(n[i], newx$MLE, logLikZ(newx), newx$nobs)
    }
    
    return(as.data.frame(o))
}

ps <- seq(1, 0.1, length.out = 9)
s1 <- oneShot(1000, 0.01, ps)
s2 <- oneShot(10000, 0.01, ps)
s1i <- itrShot(1000, 0.01, ps)
s2i <- itrShot(10000, 0.01, ps)

# plot(range(ps), range(s1$mle, s2$mle, s1i$mle, s2i$mle), type = 'n', log = 'x')
# points(ps, s1$mle, type = 'b', col = 'red', lty = 2)
# points(ps, s2$mle, type = 'b', col = 'blue', lty = 2)
# points(ps, s1i$mle, type = 'b', col = 'red', pch = 16)
# points(ps, s2i$mle, type = 'b', col = 'blue', pch = 16)

plot(range(ps), range(s1$s0, s2$s0, s1i$s0, s2i$s0), type = 'n', log = 'xy')
points(ps, s1$s0, type = 'b', col = 'red', lty = 2)
points(ps, s2$s0, type = 'b', col = 'blue', lty = 2)
points(ps, s1i$s0, type = 'b', col = 'red', pch = 16)
points(ps, s2i$s0, type = 'b', col = 'blue', pch = 16)

# plot(range(s1$n, s2$n), range(s1$mle, s2$mle), type = 'n', log = 'x')
# points(s1$n, s1$mle, type = 'b', col = 'red')
# points(s2$n, s2$mle, type = 'b', col = 'blue')

# plot(range(ps), range(s1$mle, s2$mle), type = 'n', log = 'y')
# points(ps, s1$mle, type = 'b', col = 'red')
# points(ps, s2$mle, type = 'b', col = 'blue')
# 
# plot(range(ps), range(s1$ll_z, s2$ll_z), type = 'n', log = 'y')
# abline(h = qchisq(0.975, 1))
# points(ps, s1$ll_z, type = 'b', col = 'red')
# points(ps, s2$ll_z, type = 'b', col = 'blue')

# plot(range(ps), range(s1$s0, s2$s0), type = 'n', log = 'xy')
# points(ps, s1$s0, type = 'b', col = 'red')
# points(ps, s2$s0, type = 'b', col = 'blue')

plot(range(ps), c(0.2, 1), type = 'n', log = 'xy')
points(ps, s1$s0 / max(s1$s0), type = 'b', col = 'red')
points(ps, s2$s0 / max(s2$s0), type = 'b', col = 'blue')


