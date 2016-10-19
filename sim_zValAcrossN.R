library(pika)
library(parallel)

N <- seq(10, 100, by = 15)

zn <- lapply(N, function(n) {
    out <- mclapply(1:100, mc.cores = 6, FUN = function(i) {
        print(paste(n, ': ', i, sep = ''))
        s <- sad(x = rtnegb(n, 10, 10), model = 'tpois')
        logLikZ(s, nrep = 500)$z
    })
    out <- unlist(out)
    
    return(mean = c(mean(out), quantile(out, prob = c(0.025, 0.975))))
})

zn <- data.frame(n = N, do.call(rbind, zn))

x <- sad(rtnegb(100, 10, 10), model = 'tpois', keepData = TRUE)
plot(x, ptype = 'rad')

xr <- getrfun(x)
foo <- replicate(500, sad(xr(100), 'stick')$ll)

plot(density(foo), xlim = range(foo, x$ll))
abline(v = x$ll)

foo <- logLikZ(x, nrep = 999, return.sim = TRUE)
foo$z


