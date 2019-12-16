context('expectedN functions work')


test_that('solve4expectedN and expectedN are consistent', {
    S <- 200
    N <- S:(10000 * S)
    Nhat <- expectedN(S, 'fish', 0.01)
    sol <- solve4expectedN(round(Nhat), 'fish', 0.01)

    expect_true(abs(sol$root - S) / S < 0.01)
})

test_that('dNGivenS and expectedN are consistent', {
    S <- 200
    N <- S:(10000 * S)
    p <- dNGivenS(N, S, 0.01, 'fish')

    expect_equal(sum(p * N), expectedN(S, 'fish', 0.01))
})

