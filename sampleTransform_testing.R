source('sampleTransform.R')

# x <- sad(rtnegb(500, 30, 0.5), keepData = TRUE)
bci <- read.csv('~/Dropbox/Research/data/stri/BCIS.csv', as.is = TRUE)
bci <- bci[bci$year == max(bci$year), ]
bci <- tapply(bci$count, bci$spp, sum)

x <- sad(bci, keepData = TRUE)
plot(x, ptype = 'rad', log = 'y')

J <- sum(x$data)
N <- round(exp(seq(log(0.01*J), log(0.75*J), length.out = 10)))

ll <- lapply(N, function(n) {
    xsamp <- sample.sad(x, n)
    
    return(c(llST = sum(sampTrans(xsamp, x, n, include0 = FALSE, log = TRUE)),
             llF = logLik(sad(xsamp, model = 'fish', keepData = TRUE)),
             llPLN = logLik(sad(xsamp, model = 'plnorm', keepData = TRUE)),
             llNB = logLik(sad(xsamp, model = 'tnegb', keepData = TRUE))))
    
})

ll <- do.call(rbind, ll)

matplot(N, ll - ll[, 1], type = 'l', log = 'x', ylim = c(-60, 1), lty = 1, 
        col = rainbow(ncol(ll), end = 0.7))
legend('bottomleft', legend = colnames(ll), lty = 1, col = rainbow(ncol(ll), end = 0.7))
