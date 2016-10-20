library(meteR)
library(parallel)

setwd('~/Dropbox/Research/happySAD')

N <- round(exp(seq(log(10), log(200), length = 10)) / 5) * 5

zmeteRn <- lapply(N, function(n) {
  out <- mclapply(1:100, mc.cores = 6, FUN = function(i) {
      print(paste(n, ': ', i, sep = ''))
      x <- rnbinom(100*n, mu = 8, size = 1)
      x <- sample(x[x > 0], n, replace = ifelse(sum(x > 0) >= n, FALSE, TRUE))
      s <- sad(meteESF(spp = as.character(1:n), abund = x))
      z <- logLikZ(s, nrep = 500)$z
      return(z)
  })
  out <- unlist(out)
  
  return(c(mean = mean(out), quantile(out, prob = c(0.025, 0.975))))
})

zmeteRn <- data.frame(n = N, do.call(rbind, zmeteRn))
names(zmeteRn) <- c('n', 'mean', 'ciLo', 'ciHi')

write.csv(zmeteRn, file = 'sim_zVal_meteR.csv', row.names = FALSE)
