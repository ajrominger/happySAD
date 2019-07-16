pmultinom(rep(3, 12), 12, rep(1, 12), 12)





nn <- 0:200
ll <- seq(1, 5, length.out = 20)
fooThr <- sapply(nn, dsumuptpois, lambda = ll, m = rep(100, length(ll)))
plot(nn, fooThr, ylim = c(0, max(fooThr)))
points(nn, dpois(nn, sum(ll)), col = 'red')