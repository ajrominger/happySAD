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

sad.rad <- lapply(names(sad.par), function(m) {
    lapply(1:4, function(i) {
        s <- sad(model=m, par=sad.par[[m]][[i]])
        return(sad2Rank(s, S=nspp[4]))
    })
})
names(sad.rad) <- names(sad.par)

pdf(file='ms/fig/fig_sadsUsed.pdf', width=4, height=4)
par(mfcol=c(4, 4), mar=rep(0.5, 4), oma=c(2, 4, 2, 0)+0.1, mgp=c(2, 0.75, 0))

ylim <- range(unlist(sad.rad))

for(i in 1:4) {
    mod <- names(sad.rad)[i]
    for(j in 1:4) {
        plot(sad.rad[[i]][[j]], ylim=ylim, type='l', log='y', xaxt='n', yaxt='n',
             col=cols[mod], lwd=2)
        if(i == 1) logAxis(2)
        if(j == 1) mtext(switch(mod, 'fish'='Logseries',
                                'plnorm'='PoisLogNorm',
                                'stick'='BrokenStick',
                                'tnegb'='TruncNegBin'), 
                         side=3, line=1, cex=0.8, col=cols[mod])
    }
}

mtext('Species rank', side=1, line=1, outer=TRUE)
mtext('Abundance', side=2, line=2, outer=TRUE)

dev.off()



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
pdf('ms/fig/fig_aic.pdf', width=4, height=4)
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
                                      'plnorm'='PoisLogNorm',
                                      'stick'='BrokenStick',
                                      'tnegb'='TruncNegBin'), 
                               side=3, line=1, cex=0.8, col=cols[mod])
        
        if(mod=='tnegb') mtext(pr, side=4, line=1)
    }
}

mtext('Different parameterizations', side=1, line=1.5, outer=TRUE)
mtext('Relative model support', side=2, line=1, outer=TRUE)

dev.off()


## ===============================
## plot Z-value simulation results
## ===============================

getZ <- function(stat, porz) {
    stat <- paste(porz, stat, sep='_')
    these <- sim[, 1, 1] == stat & !is.na(sim[, 'value', 1])
    stat.sim <- sim[these, 'value', ]
    samedist <- sim[these, 'fittedDist', 1] == sim[these, 'actualDist', 1]
    out <- list(right=NULL, wrong=NULL)
    out$right <- as.numeric(stat.sim[samedist, ])
    out$wrong <- as.numeric(stat.sim[!samedist, ])
    
    return(out)
}

zll <- getZ('cdfMSErel', 'z')

plot(density(log(zll$right)))
lines(density(log(zll$wrong)))


bla <- as.data.frame(sim[1:(15*4), , 1])
bla <- dcast(as.data.frame(sim[, , 1]), fittedDist + actualDist + pars + prop + nspp ~ stat)
head(bla)

sim2 <- lapply(sim.out, function(x) {
    out <- melt(x)
    colnames(out) <- c('stat', 'fittedDist', 'value', 'pars', 'actualDist', 'prop', 'nspp')
    out$pars <- pars[paste(out$actualDist, out$pars, sep='')]
    out$prop <- prop[out$prop]
    out$nspp <- nspp[out$nspp]
    out$stat <- as.character(out$stat)
    out$fittedDist <- as.character(out$fittedDist)
    
    return(dcast(out, actualDist + pars + prop + nspp + fittedDist ~ stat))
})

sim2 <- do.call(rbind, sim2)

## make summary of z score have deltaAIC 
z.sum <- sim2
for(i in seq(0, nrow(z.sum)-1, by=4)) {
    df <- z.sum[i + 1:4, ]
    refAIC <- df$aic[df$fittedDist == df$actualDist]
    z.sum$aic[i + 1:4] <- z.sum$aic[i + 1:4] - refAIC
    if(any(z.sum$aic[i+1:4] < 0)) z.sum$aic[i + 1:4] <- rep(NA, 4)
}



plot(density(z.sum[z.sum$actualDist != z.sum$fittedDist, 'aic']))

plot(z.sum[z.sum$actualDist!='tnegb' & z.sum$fittedDist=='tnegb', c('aic', 'z_radMSE')], col='red', log='y')
points(z.sum[z.sum$actualDist=='tnegb' & z.sum$fittedDist=='tnegb', c('aic', 'z_radMSE')])

plot(z.sum[z.sum$actualDist!='plnorm' & z.sum$fittedDist=='plnorm', c('aic', 'z_radMSE')])


plot(density(log(z.sum[z.sum$actualDist=='tnegb' & z.sum$fittedDist=='tnegb', 'z_radMSE'])))



x <- rplnorm(300, 0, 1.5)
plot(sort(x, TRUE), log='y')
tnegbfit <- sad(x, 'tnegb', keepData=TRUE)
plnormfit <- sad(x, 'plnorm', keepData=TRUE)
stickfit <- sad(x, 'stick', keepData=TRUE)
lines(sad2Rank(tnegbfit))
lines(sad2Rank(plnormfit), col='red')

pz <- mseZ(plnormfit, nrep=1000, return.sim=TRUE, rel=TRUE)
tz <- mseZ(tnegbfit, nrep=1000, return.sim=TRUE, rel=TRUE)
sz <- mseZ(stickfit, nrep=1000, return.sim=TRUE, rel=TRUE)

pz <- logLikZ(plnormfit, nrep=1000, return.sim=TRUE)
tz <- logLikZ(tnegbfit, nrep=1000, return.sim=TRUE)
sz <- logLikZ(stickfit, nrep=1000, return.sim=TRUE)
plot(density(log(pz$sim)))
lines(density(log(tz$sim)))
lines(density(log(sz$sim)))




bla <- z.sum[z.sum$aic < 0, 1:6]
nrow(bla)
aggregate(list(aic=bla$aic), bla[, c('actualDist', 'fittedDist', 'pars')], min)

x <- replicate(100, {
    r <- rplnorm(140, 2, 0.3)
    fit <-  fitSAD(r, c('plnorm', 'tnegb'))
    if(fit[[2]]$ll - fit[[1]]$ll < -20) browser()
})
