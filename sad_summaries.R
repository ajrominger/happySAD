library(pika) # can be installed with devtools: `devtools::install_github('ajrominger/pika')`
library(socorro) # can be installed with devtools: `devtools::install_github('ajrominger/socorro')`

#' function implementing renyi entropies for an SAD
#' @param x a vector of species abundances
#' @param alpha the entropy order defining each entropy 
#' $H_\alpha = \frac{1}{1 - \alpha} \text{log}\right(\sum_i p_i^\alpha\left)$

renyi <- function(x, alpha) {
    p <- table(x)
    p <- as.numeric(p) / sum(p)
    
    out <- 1 / (1 - alpha) * log(colSums(outer(p, alpha, '^')))
    out[alpha == 1] <- - sum(p * log(p))
    
    return(out)
}

# pars of the truncated neg binom to loop over
# we simulate from the truncated negative binomial because it has a very
# flexible shape that can look very log-series (k < 0.1), or very log-normal
# (k > 1; mu > 10 * k, roughly) and also geometric (k = 1) and Poisson (k very large)

nsim <- 500
pars <- cbind(mu = runif(nsim, 0.5, 50), k = runif(nsim, 0.01, 10))

# compute summary statistics:
# h0.25 through h16 are renyi entropies
# sd is standard deviation
# skew is skewness
# kurt is kurtosis
sumstat <- parallel::mclapply(1:nrow(pars), mc.cores = 8, FUN = function(i) {
    x <- rtnegb(100, pars[i, 1], pars[i, 2])
    
    out <- c(h = renyi(x, c(0.25, 1, 4, 16)), 
             sd = sd(x),
             skew = skew(x),
             kurt = kurt(x))
    names(out)[1:4] <- c('h0.25', 'h1', 'h4', 'h16')
    
    return(out)
})

sumstat <- do.call(rbind, sumstat)

# see which are most like each other
pdf('sad_sumstat.pdf', width = 6, height = 6)
pairs(sumstat)
dev.off()
