library(happySAD)

bb <- c(0.005, 0.05, 0.5)

nexp <- sapply(bb, function(b) {
    log(range(qNGivenS(c(0.025, 0.975), 500, b, 'fish'),
              qNGivenS(c(0.025, 0.975), 100, b, 'fish')), 2)
})
nexp <- floor(min(nexp)):ceiling(max(nexp))
nexp <- nexp[seq(1, by = 2, length.out = 4)]

nn <- 2^nexp

nexp <- nexp[!is.na(nn)]
nn <- nn[!is.na(nn)]


s <- 1:max(nn)

allB <- lapply(bb, function(b) {
    M <- sapply(nn, function(n) {
        dSGivenN(s, n, b, 'fish')
    })
    M <- t(M)
    goodM <- colSums(M > .Machine$double.eps^0.5) > 0
    M <- M[, goodM]
    s <- s[goodM]


    P <- parallel::mclapply(1:length(s), mc.cores = 10, FUN = function(i) {
        cat(b, ': ', i, '\n', sep = '')
        dintModeGivenS(s[i], pars = b, mod = dfish)
    })

    P <- do.call(rbind, P)

    return(t(M %*% P))
})

