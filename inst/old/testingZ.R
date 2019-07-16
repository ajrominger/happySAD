DAT <- replicate(100, rnorm(200))

x <- t(apply(DAT, 2, function(X) {
  dat <- X
#   dat <- scale(dat)[, 1]
  sw.p <- shapiro.test(dat)$p.value
  ks.p <- ks.test(dat, pnorm)$p.value
  
  obs.lik <- sum(dnorm(dat, log=TRUE))
  sim.lik <- replicate(299, {
    sum(dnorm(rnorm(200), log=TRUE))
  })
  sim.lik <- c(obs.lik, sim.lik)
  
  lz.p <- sum(sim.lik < obs.lik)/length(sim.lik)
  
  c(sw.p, ks.p, lz.p)
}))

plot(x[, c(1, 3)], xlim=range(x, na.rm=TRUE), ylim=range(x, na.rm=TRUE)); abline(0, 1)
