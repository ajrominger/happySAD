context('distribution functions work')

# ----
# dsumuptpois

test_that('when upper truncation is large the `dsumuptpois` function is equivilant to un-truncated', {
    nn <- 0:200
    ll <- seq(1, 5, length.out = 20)
    upt <- sapply(nn, dsumuptpois, lambda = ll, m = rep(10000, length(ll)))
    noup <- dpois(nn, sum(ll))

    expect_true(mean(abs(upt - noup)) < 1e-03)
    expect_true(abs(sum(upt) - sum(noup)) < 1e-03)
})


# ----
# pmultinom

test_that('pmultinom output matches published value', {
    levin <- 0.8367
    hyp <- pmultinom(rep(3, 12), 12, rep(1, 12), 12)
    expect_equal(round(hyp * 1000), round(levin * 1000))
})


# ----
# dintModeGivenS


# b <- 10
# S <- 1000
# noct <- 13
# nn <- 1:(2^(noct) - 1)
# octs <- floor(log(nn, 2))
#
# pp <- dfish(nn, b)
# octProbs <- tapply(pp, octs, sum)
#
# x <- proc.time()
# thr <- dintModeGivenS(S, b, dfish)
# proc.time() - x
#
# x <- proc.time()
# rr <- rmultinom(10000, S, octProbs)
# sim <- sapply(1:13, function(i) {
#     o <- apply(rr, 2, function(x) all(x[i] >= x[-i]))
#     mean(o)
# })
# proc.time() - x
#
# plot(thr, ylim = 0:1)
# points(sim, col = 'red')
#
#

#
# N <- 500
# ss <- 10:150
# b <- 0.01
# thr <- dSGivenN(ss, N, b, 'fish')
#
#
#
# nrep <- 10000
# x <- rfish(max(ss) * nrep, b)
# sim <- parallel::mclapply(ss, mc.cores = 3, FUN = function(s) {
#     # x <- rfish(s * nrep, b)
#     # xx <- split(x, rep(1:nrep, each = s))
#     # sapply(xx, sum)
#     #
#     replicate(nrep, {
#         sum(sample(x, s, replace = TRUE))
#     })
# })
#
# sim <- data.frame(s = rep(ss, each = nrep), N = unlist(sim))
# sim <- as.data.frame(simpECDF(sim$s[sim$N == N]))
# names(sim) <- c('x', 'psim')
#
# thr <- data.frame(x = ss, pthr = cumsum(thr))
#
# plot(thr, log = 'xy')
# points(sim, col = 'red')
#
# comp <- merge(thr, sim, by = 'x')
# sum((comp$pthr - comp$psim)^2)
