## =====================================================================
## script to loop over several values for number of species, simulate a
## truncated negative binomial and fit the METE SAD to it, calculate the
## z value and see how that changes with number of species
## =====================================================================

library(meteR)
library(parallel)

setwd('~/Dropbox/Research/happySAD')

## numbers of species
N <- round(exp(seq(log(10), log(200), length = 10)) / 5) * 5

## loop over them
zmeteRn <- lapply(N, function(n) {
  out <- mclapply(1:100, mc.cores = 6, FUN = function(i) {
      print(paste(n, ': ', i, sep = ''))
      
      ## simulate truncated neg binom
      x <- rnbinom(100*n, mu = 8, size = 1)
      x <- sample(x[x > 0], n, replace = ifelse(sum(x > 0) >= n, FALSE, TRUE))
      
      ## fit sad
      s <- sad(meteESF(spp = as.character(1:n), abund = x))
      
      ## calculate Z
      z <- logLikZ(s, nrep = 500)$z
      return(z)
  })
  out <- unlist(out)
  
  return(c(mean = mean(out), quantile(out, prob = c(0.025, 0.975))))
})


## clean-up output and write out
zmeteRn <- data.frame(n = N, do.call(rbind, zmeteRn))
names(zmeteRn) <- c('n', 'mean', 'ciLo', 'ciHi')
write.csv(zmeteRn, file = 'sim_zVal_meteR.csv', row.names = FALSE)


## plot it
pdf('figs/fig_zValN_meteR.pdf', width = 5, height = 5)
par(mar = c(3, 3, 0, 0) + 0.5, mgp = c(2, 0.75, 0))
plot(zmeteRn[, 1:2], ylim = range(zmeteRn[, -1]), 
     xlab = 'Number of species', ylab = expression(z^2~'value'),
     pch = 21, bg = 'white', cex = 1.2,
     panel.first = arrows(x0 = zmeteRn[, 1], y0 = zmeteRn$ciLo, y1 = zmeteRn$ciHi,
                          code = 3, angle = 90, length = 0.05))
dev.off()
