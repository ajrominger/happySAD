library(pika)
bb <- 1
nprop <- seq(1, 1, length.out = 1)
nsp <- 8

foo <- lapply(bb, function(b) {
    x <- rpois(100000, b)
    # print(any(!is.finite(x)))
    x[!is.finite(x)] <- max(x[is.finite(x)])
    sadList <- split(x, ceiling(seq_along(x) / nsp))
    
    samps <- sapply(nprop, function(np) {
        intMode <- sapply(sadList, function(thisOne) {
            # n <- round(np * sum(thisOne))
            # thisOne <- sample.sad(sad(thisOne, keepData = TRUE), n, replace = FALSE)
            # octCount <- table(floor(log(thisOne, 2)))
            
            octCount <- table(thisOne)
            
            return(any(octCount[1] < octCount[2]))
        })

        return(mean(intMode, na.rm = TRUE))
    })

    return(samps)
})
foo

dintModeGivenS <- function(S, b) {
    browser()
    nn <- 1:2
    octProbs <- dpois(nn, b)
    
    y <- 1:S
    
    o <- pbinom(floor((y - 1) / 2), y, octProbs[1] / sum(octProbs[c(1, 2)])) *
        dbinom(y, S, sum(octProbs[c(1, 2)]))
    
    return(sum(o))
}

foo
dintModeGivenS(nsp, bb)

dintModeGivenS <- function(S, b) {
    nn <- 1:(2^13 - 1)
    octs <- floor(log(nn, 2))
    probs <- rfish(nn, b)
    octProbs <- tapply(probs, octs, sum)
    octProbs <- c(octProbs, '13' = 0)
    octProbs <- octProbs / sum(octProbs)
    p0 <- octProbs[1]
    
    dS0i <- function(S0i, i) {
        dbinom(S0i, S, p0 + octProbs[i])
    }
    
    db0lessi <- function(S0i, i) {
        pbinom(floor(S0i / 2), S0i, p0 / (p0 + octProbs[i]))
    }
    
    y <- 1:S
    o <- sapply(2:length(octProbs), function(i) {
        db0lessi(y, i) * dS0i(y, i)
    })
    
    return(sum(o))
}

dintModeGivenS(nsp, bb)


dintModeGivenS <- function(S, b) {
    # browser()
    nn <- 1:(2^13 - 1)
    octs <- floor(log(nn, 2))
    probs <- dfish(nn, b)
    octProbs <- tapply(probs, octs, sum)
    octProbs <- c(octProbs, '13' = 0)
    octProbs <- octProbs / sum(octProbs)

    y <- 1:S

    # o <- sapply(2:length(octProbs), function(i) {
    #     pbinom(floor(y / 2), y, octProbs[1] / sum(octProbs[c(1, i)])) *
    #         dbinom(y, S, sum(octProbs[c(1, i)]))
    # })
    
    o <- pbinom(floor(y / 2), y, octProbs[1] / sum(octProbs[c(1, 2)])) *
            dbinom(y, S, sum(octProbs[c(1, 2)]))

    return(sum(o))
}







