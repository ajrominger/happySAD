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
            out <- function(n) get(sprintf('r%s', f))(n, sad.par[[f]][[p]][1], sad.par[[f]][[p]][2])
        } else {
            out <- function(n) get(sprintf('r%s', f))(n, sad.par[[f]][p])
        }
        
        attr(out, 'model') <- f
        out
    })
})

names(sad.rfun) <- names(sad.par)


## ==============================
## function to run one simulation
## ==============================

simSAD <- function(funs, nspp, prop) {
    rapply(funs, how='replace', f=function(f) {
        ## simulate data
        #dat <- try(sample.sad(f(nspp), prob=prop))
        dat <- sample.sad(f(nspp), prob=prop)
        #if(class(dat) == 'try-error') {
        #return(matrix(1, nrow=3, ncol=4))
        #} else {
        ## fit SAD models and extract AIC
        #fit <- try(fitSAD(dat, c('fish', 'plnorm', 'stick', 'tnegb'), keepData=FALSE))
        fit <- fitSAD(dat, c('fish', 'plnorm', 'stick', 'tnegb'), keepData=FALSE)
        #if(class(fit) == 'try-error') {
        #return(matrix(2, nrow=3, ncol=4))
        #} else {
        aic <- sapply(fit, AIC)
        aicWin <- aic - min(aic) <= 2
        aicWin <- aicWin/sum(aicWin)
        daic <- aic - aic[attr(f, 'model')]
        
        return(rbind(aic=aic, aicWin=aicWin, daic=daic))
        #}
        #}
    }) 
}

## ==============
## run simulation
## ==============

nsim <- 500
if(Sys.info()['nodename'] == 'voldemort') {
    mc <- 12
} else {
    mc <- 6
}

sim.out.aic <- mclapply(1:nsim, mc.cores=mc, FUN=function(n) {
    cat(n, '\n')
    lapply(nspp, function(ns) {
        lapply(prop, function(p) {
            simSAD(sad.rfun, ns, p)
        })
    })
})


## save simulation output
save(sim.out.aic, nspp, prop, sad.par, sad.rfun, file=sprintf('sim_out_aic_%s.RData', Sys.info()['nodename']))
