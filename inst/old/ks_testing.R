N <- c(10000, seq(100, 2000, length.out = 10))

ksZ <- lapply(N, function(n) {
    nullStat <- replicate(100, {
        x2 <- rgamma(n, 1, 0.7)
        ks.test(x2, 'pgamma', 1, 0.7)$statistic
    })
    
    out <- replicate(50, {
        x1 <- rlnorm(n, 1.5, 0.6)
        stat <- ks.test(x1, 'pgamma', 1, 0.7)$statistic
        
        return((stat - mean(nullStat)) / sd(nullStat))
    })
    
    c(mean(out), quantile(out, c(0.025, 0.975)))
})

ksZ <- do.call(rbind, ksZ)

ks <- lapply(N, function(n) {
    out <- replicate(50, {
        x1 <- rnbinom(n, size = 5, mu = 3)
        ks.test(x1, 'rnbinom', size = 0.001, mu = 4)$statistic
    })
    
    c(mean(out), quantile(out, c(0.025, 0.975)))
})

ks <- do.call(rbind, ks)


par(mfrow = 1:2, mar = c(2, 2, 0, 0) + 0.5)
curve(dlnorm(x, 1.5, 0.6), from = 0, to = 10)
curve(dgamma(x, 1, 0.7), add = TRUE, col = 'red')
plot(N[-1], ksZ[-1, 1], ylim = range(ksZ), pch = 16,
     panel.first = {
         rect(xleft = par('usr')[1], xright = par('usr')[2], 
              ybottom = ksZ[1, 2], ytop = ksZ[1, 3], 
              col = 'gray', border = NA)
         abline(h = ksZ[1, 1])
         segments(x0 = N[-1], y0 = ksZ[-1, 2], y1 = ksZ[-1, 3])
     }
)

par(mfrow = 1:2, mar = c(2, 2, 0, 0) + 0.5)
curve(dlnorm(x, 1.5, 0.6), from = 0, to = 10)
curve(dgamma(x, 1, 0.7), add = TRUE, col = 'red')
plot(N[-1], ks[-1, 1], ylim = range(ks), pch = 16,
     panel.first = {
         rect(xleft = par('usr')[1], xright = par('usr')[2], 
              ybottom = ks[1, 2], ytop = ks[1, 3], 
              col = 'gray', border = NA)
         abline(h = ks[1, 1])
         segments(x0 = N[-1], y0 = ks[-1, 2], y1 = ks[-1, 3])
     }
)
