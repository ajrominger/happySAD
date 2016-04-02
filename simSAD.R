## ======================================================================
## simulation that generates random SAD data under different conditions:
##   1) different numbers of species
##   2) different subsampling sizes
##   3) different models
##   4) different parameter values
## and then tests:
##   1) AIC methods
##   2) likelihood goodness of fit methods
##   3) mse goodness of fit methods on both rank and cumulative functions
## ======================================================================

setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')


## =================
## simulation set up
## =================

## number of species
nspp <- 10^seq(1.75, 2.75, length=5) - 10^seq(1.75, 2.75, length=5) %% 50

## proportion sampled
prop <- 10^seq(-0.6, 0, length=5) - 10^seq(-0.6, 0, length=5) %% 0.05

## SAD parameters
sad.par <- list(fish=10^seq(0, -2.5, length=4),
                plnorm=list(c(0, 1), c(3, 1), c(1, 2), c(2.5, 2.5)),
                stick=seq(0.5, 1/100, length=4),
                tnegb=list(c(1, 1), c(8, 1), c(20, 1.5), c(100, 1.5)),
                tpois=seq(2, 15, length=4))


## ======================
## save simulation output
## ======================
save(nspp, prop, sad.par, file='sim_out.RData')
