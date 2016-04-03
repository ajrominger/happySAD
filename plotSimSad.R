## ======================================================================
## plotting SAD simulation
## ======================================================================


setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
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
        
        # legend('topright', legend=mean(rfun(100, sad.par[[i]][j])))
    }
    
    mtext(names(sad.par[i]), side=4, line=1, xpd=NA)
}

mtext('Abundance', side=1, line=3, outer=TRUE)
mtext('Probability', side=2, line=2, outer=TRUE)
