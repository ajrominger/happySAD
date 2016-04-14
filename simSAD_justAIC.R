## ======================================================================
## simulation that generates random SAD data under different conditions:
##   1) different numbers of species
##   2) different subsampling sizes
##   3) different models
##   4) different parameter values
## and then tests:
##   1) AIC methods
## ======================================================================

setwd('~/Dropbox/Research/happySAD')
# devtools::load_all('../pika')
install.packages('~/Dropbox/Research/pika', repos=NULL, type='source')
# devtools::install_github('ajrominger/pika')
library(pika)
library(parallel)


## =================
## simulation set up
## =================

## number of species
nspp <- 10^seq(1.95, 2.7, length=4)
nspp <- nspp - nspp %% 20

## proportion sampled
prop <- 10^seq(-0.4, 0, length=4)
prop <- prop - prop %% 0.05

## SAD parameters
sad.par <- list(fish=10^seq(-1.25, -2, length=4),
                plnorm=list(c(2, 0.5), c(3, 0.5), c(0, 1.5), c(1, 1.5)),
                stick=10^seq(-0.7, -1.7, length=4),
                tnegb=list(c(20, 4), c(8, 1.5), c(8, 0.5), c(100, 3))
                # tpois=seq(5, 50, length=4)
)

sad.rfun <- lapply(names(sad.par), function(f) {
    lapply(1:4, function(p) {
        if(f %in% c('plnorm', 'tnegb')) {
            return(function(n) get(sprintf('r%s', f))(n, sad.par[[f]][[p]][1], sad.par[[f]][[p]][2]))
        } else {
            return(function(n) get(sprintf('r%s', f))(n, sad.par[[f]][p]))
        }
    })
})

names(sad.rfun) <- names(sad.par)


## ==============================
## function to run one simulation
## ==============================

simSAD <- function(rfuns, nspp, prop) {
    rapply(rfuns, how='replace', f=function(f) {
        ## simulate data
        dat <- sample.sad(f(nspp), prob=prop)
        
        ## fit SAD models and extract AIC
        fit <- fitSAD(dat, c('fish', 'plnorm', 'stick', 'tnegb'), keepData=FALSE)
        return(sapply(fit, AIC))
    }) 
}

## ==============
## run simulation
## ==============

nsim <- 500
sim.out.aic <- mclapply(1:nsim, mc.cores=16, FUN=function(n) {
    print(n)
    lapply(nspp, function(ns) {
        lapply(prop, function(p) {
            simSAD(sad.rfun, ns, p)
        })
    })
})

## save simulation output
save(sim.out.aic, nspp, prop, sad.par, sad.rfun, file='sim_out_aic.RData')
