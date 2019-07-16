library(pika)
library(parallel)

N <- round(exp(seq(log(10), log(200), length = 10)) / 5) * 5

zn <- lapply(N, function(n) {
    out <- mclapply(1:100, mc.cores = 6, FUN = function(i) {
        print(paste(n, ': ', i, sep = ''))
        s <- sad(x = rtnegb(n, 10, 20), model = 'tpois', keepData = TRUE)
        z <- logLikZ(s, nrep = 500)$z
        e <- mseZ(s, nrep = 500, type = 'rank', relative = FALSE, log = FALSE)$z
        return(c(z = z, e = e))
    })
    out <- do.call(rbind, out)
    
    return(c(z = c(mean = mean(out[, 1]), quantile(out[, 1], prob = c(0.025, 0.975))),
             e = c(mean = mean(out[, 2]), quantile(out[, 2], prob = c(0.025, 0.975)))))
})

zn <- data.frame(n = N, do.call(rbind, zn))
names(zn) <- c('n', 'z.mean', 'z.ciLo', 'z.ciHi', 'e.mean', 'e.ciLo', 'e.ciHi')

plot(zn$n, zn$e.mean)


## now look at when fitted model and generating model are the same

zn <- lapply(N, function(n) {
    out <- mclapply(1:100, mc.cores = 6, FUN = function(i) {
        print(paste(n, ': ', i, sep = ''))
        s <- sad(x = rtnegb(n, 10, 20), model = 'tnegb', keepData = TRUE)
        z <- logLikZ(s, nrep = 500)$z
        return(z)
    })
    out <- unlist(out)
    
    return(c(mean = mean(out), quantile(out, prob = c(0.025, 0.975))))
})

zn <- data.frame(n = N, do.call(rbind, zn))
names(zn) <- c('n', 'mean', 'ciLo', 'ciHi')

plot(zn[, 1:2])

