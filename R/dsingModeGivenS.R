#' @title Probability of a single mode
#'
#' @description For a given model and number of species, returns the probability of a
#' single mode in the Preston plot
#'
#' @details A mode is defined by the cuttoff \code{modeCutOff} and so a mode must be greater
#' than its neighbors by at least that value
#'
#' @param x object to Prestonize
#'
#' @return a numeric giving the desired probability
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

dsingModeGivenS <- function(S, mod, pars, modeCutOff = 1, B = 1e+05) {
    noct <- 13
    x <- 1:noct
    nn <- 1:(2^(noct) - 1)
    octs <- floor(log(nn, 2))
    pp <- do.call(mod, c(list(nn), as.list(pars)))
    octProbs <- tapply(pp, octs, sum)

    r <- rmultinom(B, S, octProbs)
    # plot(r, type = 'b')
    # abline(h = seq(0, 80, by = modeCutOff))

    # browser()
    # o <- parallel::mclapply(1:ncol(r), mc.cores = 3, FUN = function(i) {
    o <- lapply(1:ncol(r), function(i) {
        .allLessMax(r[, i], modeCutOff)
    })

    return(mean(unlist(o)))
}

.allLessMax <- function(x, cutOff) {
    imax <- which.max(x)

    if(cutOff < 1) cutOff <- max(x) * cutOff
    cat(cutOff, '\n')

    maxsplit <- list(x[1:imax], x[imax:length(x)])

    return(all(diff(maxsplit[[1]]) >= 1 - cutOff) &
               all(diff(maxsplit[[2]]) <= -1 + cutOff))
}


