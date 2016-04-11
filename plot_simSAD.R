## ======================================================================
## plotting SAD simulation
## ======================================================================


setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(reshape2)
source('~/R_functions/logAxis.R')

## read in simulation
load('sim_out.RData')

## plot theoretical SADs
par(mfrow=c(5, 4), mar=rep(0.5, 4), oma=c(4, 3, 1, 3)+0.1)

for(i in 1:length(sad.par)) {
    if(names(sad.par[i]) %in% c('plnorm', 'tnegb')) {
        dfun <- function(x, par) get(sprintf('d%s', names(sad.par[i])))(x, par[[1]][1], par[[1]][2])
        rfun <- function(n, par) get(sprintf('r%s', names(sad.par[i])))(n, par[[1]][1], par[[1]][2])
    } else {
        dfun <- function(x, par) get(sprintf('d%s', names(sad.par[i])))(x, par)
        rfun <- function(n, par) get(sprintf('r%s', names(sad.par[i])))(n, par)
    }
    
    ylim <- range(sapply(sad.par[[i]], function(p) {
        if(length(p) > 1) p <- list(p)
        dfun(c(1, 50), p)
    }))
    
    
    for(j in 1:4) {
        plot(dfun(1:50, sad.par[[i]][j]), log='xy', type='l', xaxt='n', yaxt='s', lwd=3, ylim=ylim)
        if(i == 5) logAxis(1)
        
        legend('topright', legend=mean(rfun(100, sad.par[[i]][j])))
    }
    
    mtext(names(sad.par[i]), side=4, line=1, xpd=NA)
}

mtext('Abundance', side=1, line=3, outer=TRUE)
mtext('Probability', side=2, line=2, outer=TRUE)


par(mfrow=c(length(sad.par), 4), mar=c(1.5, 1.5, 0, 0)+0.1, mgp=c(2, 0.4, 0), tcl=-0.4, oma=c(4, 4, 0, 0))
for(i in 1:length(sad.par)) {
    for(j in 1:4) {
        r <- sort(sample.sad(sad.rfun[[i]][[j]](nspp[4]), prob=prop[4]), TRUE)
        plot(r, xlab='', ylab='', log='y', type='l', ylim=c(1, max(r)))
        for(k in 1:3) lines(sort(sample.sad(sad.rfun[[i]][[j]](nspp[4]), prob=prop[k]), TRUE))
    }
}


## clean simulation results

sim.out <- sim.out[sapply(sim.out, class) == 'list']

pars <- unlist(lapply(sad.par, function(x) {
    if(class(x) == 'numeric') {
        return(as.character(round(x, 3)))
    } else {
        return(sapply(x, function(X) paste(round(X, 3), collapse=', ')))
    }
}))


## combine sim results into a 3D array

sim <- lapply(sim.out, function(x) {
    out <- melt(x)
    colnames(out) <- c('stat', 'fittedDist', 'value', 'pars', 'actualDist', 'prop', 'nspp')
    out$pars <- pars[paste(out$actualDist, out$pars, sep='')]
    out$prop <- prop[out$prop]
    out$nspp <- nspp[out$nspp]
    out$stat <- as.character(out$stat)
    out$fittedDist <- as.character(out$fittedDist)
    
    return(out)
})

sim <- array(unlist(sim), dim=c(dim(sim[[1]]), length(sim)),
             dimnames=list(rownames(sim[[1]]), colnames(sim[[1]]), 1:length(sim)))


