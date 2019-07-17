context('distribution functions work')

# ----
# dsumuptpois

test_that('when upper truncation is large the `dsumuptpois` function is equivilant to un-truncated', {
    nn <- 0:200
    ll <- seq(1, 5, length.out = 20)
    upt <- sapply(nn, dsumuptpois, lambda = ll, m = rep(1000, length(ll)))
    noup <- dpois(nn, sum(ll))

    expect_true(max(abs(upt - noup)) < .Machine$double.eps^0.25)
    expect_true(abs(sum(upt) - sum(noup)) <= .Machine$double.eps^0.75)
})


# ----
# pmultinom

test_that('pmultinom output matches published value', {
    levin <- 0.8367
    hyp <- pmultinom(rep(3, 12), 12, rep(1, 12), 12)
    expect_equal(floor(hyp * 1000), floor(levin * 1000))
})
