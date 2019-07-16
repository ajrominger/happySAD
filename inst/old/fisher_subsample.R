library(pika)

x <- sad(rtnegb(500, 10, 0.2), model = 'fish', keepData = TRUE)
# x <- sad(rfish(500, 0.01), model = 'fish', keepData = TRUE)
plot(x, ptype = 'rad', log = 'y')

N <- floor(exp(seq(log(100), log(sum(x$data)), length.out = 11)))[-11]

sadsub <- lapply(N, function(n) {
    s <- sad(sample.sad(x, n), model = 'fish', keepData = TRUE)
    logLikZ(s)$z
})

plot(N, unlist(sadsub), log = 'xy')
abline(h = qchisq(0.95, 1), v = c(500, 1000))

s <- rowSums(outer(1:x$nobs, sample(rep(1:x$nobs, x$data), 500), '=='))
sum(dpois(s, lambda = x$data, log = TRUE))
logLik(sad(s[s > 0], model = 'fish', keepData = TRUE))
