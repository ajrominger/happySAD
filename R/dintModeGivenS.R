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

    return(out)
}

