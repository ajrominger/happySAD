b <- 0.1

qNGivenS(c(0.025, 0.975), 500, b, 'fish')
qNGivenS(c(0.025, 0.975), 100, b, 'fish')

nn <- 2^(6:12)



s <- 1:max(nn)
M <- sapply(nn, function(n) {
    dSGivenN(s, n, b, 'fish')
})
M <- t(M)
goodM <- colSums(M > .Machine$double.eps) > 0
M <- M[, goodM]
s <- s[goodM]

P <- parallel::mclapply(1:length(s), mc.cores = 10, FUN = function(i) {
    print(i)
    dintModeGivenS(s[i], pars = b, mod = dfish)
})

P <- do.call(rbind, P)

M %*% P
