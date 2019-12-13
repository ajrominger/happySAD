context('local maxima function works')

test_that('local maxima finder works', {
    x <- c(4, 4, 1, 2, 3, 1)
    expect_equal(localMaxima(x, by = 2, xmin = 0), c(1, 4))
})

