library(pika)
library(parallel)

setwd('~/Dropbox/Research/happySAD')

muRange <- seq(2, 10, length.out = 6)
sRange <- seq(50, 500, length.out = 6)
kRange <- seq(0.01, 1, length.out = 5)
thisTrue <- median(1:length(kRange))
k <- kRange[thisTrue]

pars <- as.matrix(expand.grid(s = sRange, mu = muRange))

allZ <- mclapply(1:nrow(pars), mc.cores = 6,
                 function(i) {
    x <- rtnegb(pars[i, 1], pars[i, 2], k)
    
    unlist(lapply(kRange, function(K) {
        thisSAD <- sad(model = 'tnegb', par = c(pars[i, 2], K))
        thisSAD$data <- x
        thisSAD$ll <- sum(dtnegb(x, pars[i, 2], K, log = TRUE))
        thisSAD$nobs <- pars[i, 1]
        
        z <- logLikZ(thisSAD, nrep = 500, return.sim = FALSE)$z
        
        return(z)
    }))
})

allZ <- do.call(rbind, allZ)
colnames(allZ) <- paste0('k', 1:length(kRange))
out <- data.frame(pars, trueK = paste0('k', thisTrue), allZ)

write.csv(out, file = 'sim_NS-ramp_zVal.csv', row.names = FALSE)
