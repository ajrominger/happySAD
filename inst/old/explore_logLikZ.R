r <- rtpois(500, 40)
rfit <- fitSAD(r, c('plnorm', 'stick', 'tnegb', 'tpois'), keepData=TRUE)
rfit[[2]]$ll
sum(dstick(rstick(500, 1/40), 1/40, log=TRUE))

plnZ <- logLikZ(rfit[[1]], return.sim=TRUE)
plot(density(plnZ$sim))
abline(v=plnZ$obs)

stickZ <- logLikZ(rfit[[2]], return.sim=TRUE)
plot(density(stickZ$sim))
abline(v=stickZ$obs)

tnbZ <- logLikZ(rfit[[3]], return.sim=TRUE)
plot(density(tnbZ$sim))
abline(v=tnbZ$z)

tpoZ <- logLikZ(rfit[[4]], return.sim=TRUE)
plot(density(tpoZ$sim))
abline(v=tpoZ$obs)

plot(.ecdf(r))
curve(getpfun(rfit[[2]])(x), col='red', add=TRUE)
