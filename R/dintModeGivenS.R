dintModeGivenS <- function(S, pars, mod) {
    noct <- 13
    x <- 1:noct
    nn <- 1:(2^(noct) - 1)
    octs <- floor(log(nn, 2))
    pp <- do.call(mod, c(list(nn), as.list(pars)))
    octProbs <- tapply(pp, octs, sum)


    possModes <- ceiling(S / noct):S
    ii <- range(which(dbinom(possModes, S, max(octProbs)) > .Machine$double.eps),
                which(dbinom(possModes, S, min(octProbs)) > .Machine$double.eps))

    # limit the number of things we're integrating over
    likModes <- possModes[ii[1]:ii[2]]
    goodx <- dbinom(min(likModes), S, octProbs) > .Machine$double.eps

    # browser()
    o <- sapply(x[goodx], function(i) {
        allS <- sapply(likModes, function(thisS) {
            probSinBi <- dbinom(thisS, S, octProbs[i])
            probRestLessS <- pmultinom(rep(thisS, noct - 1), S - thisS, octProbs[-i])
            return(probSinBi * probRestLessS)
        })

        return(sum(allS))
    })

    out <- numeric(noct)
    out[goodx] <- o

    # pmultinom(rep(thisS, noct), S - thisS, octProbs[-i])
    # mean(apply(rmultinom(10000, S - thisS, octProbs[-i]), 2, function(x) all(x <= thisS)))
    #
    # rr <- rmultinom(10000, S, octProbs)
    # doo <- sapply(1:13, function(i) {
    #     o <- apply(rr, 2, function(x) all(x[i] >= x[-i]))
    #     mean(o)
    # })
    # plot(doo, ylim = c(0, 1))
    # points(foo, col = 'red')
    return(out)
}


b <- 0.01
S <- 50
noct <- 13
nn <- 1:(2^(noct) - 1)
octs <- floor(log(nn, 2))

pp <- dfish(nn, b)
octProbs <- tapply(pp, octs, sum)

x <- proc.time()
thr <- dintModeGivenS(S, b, dfish)
proc.time() - x

x <- proc.time()
rr <- rmultinom(10000, S, octProbs)
sim <- sapply(1:13, function(i) {
    o <- apply(rr, 2, function(x) all(x[i] >= x[-i]))
    mean(o)
})
proc.time() - x

plot(thr, ylim = 0:1)
points(sim, col = 'red')


