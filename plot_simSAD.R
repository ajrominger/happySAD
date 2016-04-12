## =======================
## plotting SAD simulation
## =======================

## set-up
setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(reshape2)
source('~/R_functions/logAxis.R')

## read in simulation
load('sim_out.RData')

## standarize colors for models
cols <- hsv(c(0.52, 0.02, 0.12, 0.65), c(0.8, 0.7, 0.8, 0.7), c(0.8, 0.75, 0.9, 0.6))
names(cols) <- c('fish', 'plnorm', 'stick', 'tnegb')


## =====================
## plot theoretical SADs
## =====================

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


## ========================
## clean simulation results
## ========================

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


## ===========================
## plot AIC simulation results
## ===========================

## extract and summarize AIC

sim.aic <- sim[seq(1, dim(sim)[1], by=15), -1, ]

aic.sum <- lapply(seq(0, dim(sim.aic)[1]-1, by=4), function(i) {
    s <- sim.aic[i+(1:4), 2, ]
    
    awin <- t(apply(s, 2, function(x) {
        x <- as.numeric(x)
        wins <- x - min(x) < 2
        wins/sum(wins)
    }))
    
    ainfo <- as.data.frame(t(sim.aic[i+1, c('actualDist', 'pars', 'prop'), ]))
    
    out <- cbind(ainfo, awin)
    colnames(out) <- c(colnames(out)[1:3], c('fishWin', 'plnormWin', 'stickWin', 'tnegbWin'))
    
    return(out)
})

aic.sum <- do.call(rbind, aic.sum)
aic.sum <- aggregate(aic.sum[, c('fishWin', 'plnormWin', 'stickWin', 'tnegbWin')], aic.sum[, c('actualDist', 'pars', 'prop')], sum)


## plot AIC
par(mfcol=c(4, 4), mar=rep(0.5, 4), oma=c(3, 3, 2, 2)+0.1, mgp=c(1, 0.1, 0))

for(mod in as.character(unique(aic.sum$actualDist))) {
    for(pr in as.character(unique(aic.sum$prop))) {
        if(pr == '1') names.arg <- 1:4
        else names.arg <- rep(NA, 4)
        
        barplot(t(as.matrix(aic.sum[aic.sum$actualDist==mod & aic.sum$prop==pr, 4:7])),
                col=(cols), axes=FALSE, names.arg=names.arg, 
                space=0, xaxs='i')
        box()
        
        if(pr == '0.35') mtext(switch(mod, 'fish'='Logseries',
                                      'plnorm'='Pois LogNorm',
                                      'stick'='Broken Stick',
                                      'tnegb'='Trunc NegBin'), 
                               side=3, line=1, cex=0.8, col=cols[mod])
        
        if(mod=='tnegb') mtext(pr, side=4, line=1)
    }
}

mtext('Different parameterizations', side=1, line=1.5, outer=TRUE)
mtext('Proportion model support', side=2, line=1, outer=TRUE)

