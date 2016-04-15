## set up
setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(reshape2)

## read in files
files <- list.files('.', pattern='sim_out_aic')
sim.aic <- c()
for(f in files) {
    load(f)
    sim.aic <- c(sim.aic, sim.out.aic)
}

## clean simulation
aic.sim <- melt(sim.aic)
colnames(aic.sim) <- c('stat', 'fittedDist', 'value', 'pars', 'actualDist', 'prop', 'nspp', 'rep')
aic.sim <- dcast(aic.sim, actualDist + pars + fittedDist + nspp + prop + rep ~ stat, value.var='value')


## plot aic wins
aicWin <- with(aic.sim, aggregate(list(aicWin=aicWin), aic.sim[, c('actualDist', 'pars', 'prop', 'fittedDist')], sum))
aicWin <- dcast(aicWin, actualDist + pars + prop ~ fittedDist, value.var='aicWin')
dim(aicWin)

pdf('ms/fig/fig_aic.pdf', width=4, height=4)
par(mfcol=c(4, 4), mar=rep(0.5, 4), oma=c(3, 3, 2, 2)+0.1, mgp=c(1, 0.1, 0))

for(mod in as.character(unique(aicWin$actualDist))) {
    for(pr in 1:4) {
        if(pr == 4) names.arg <- 1:4
        else names.arg <- rep(NA, 4)
        
        barplot(t(as.matrix(aicWin[aicWin$actualDist==mod & aicWin$prop==pr, 4:7])),
                col=(cols), axes=FALSE, names.arg=names.arg, 
                space=0, xaxs='i')
        box()
        
        if(pr == 1) mtext(switch(mod, 'fish'='Logseries',
                                 'plnorm'='PoisLogNorm',
                                 'stick'='BrokenStick',
                                 'tnegb'='TruncNegBin'), 
                          side=3, line=1, cex=0.8, col=cols[mod])
        
        if(mod=='tnegb') mtext(prop[pr], side=4, line=1)
    }
}

mtext('Different parameterizations', side=1, line=1.5, outer=TRUE)
mtext('Relative model support', side=2, line=1, outer=TRUE)

dev.off()



## plot deltaAIC
mods <- as.character(unique(aicWin$actualDist))

pdf('ms/fig/fig_daic.pdf', width=4, height=4)
par(mfcol=c(4, 4), mar=c(0.1, 0.5, 0.1, 0.5), oma=c(3, 3, 2, 2)+0.1, mgp=c(1, 0.75, 0))

for(a in mods) {
    for(pr in 1:4) {
        if(a == 'tnegb') {
            plot(density(aic.sim[aic.sim$actualDist == a & aic.sim$prop==pr & aic.sim$fittedDist == 'plnorm', 'daic']), 
                 col=cols['plnorm'], yaxt='n', xaxt='n', main='', xlim=c(-10, 30), lwd=3, 
                 panel.first=rect(par('usr')[1], par('usr')[3], 2, par('usr')[4], col='gray', border=NA))
            mtext(prop[pr], side=4, line=1)
            aat <- c(-10, 10, 30)
        } else {
            if(a %in% c('fish', 'stick')) {
                xlim <- c(-1, 2.5)
                aat <- c(0, 1, 2)
            } else {
                xlim <- c(-10, 30)
                aat <- c(-10, 10, 30)
            }
            plot(density(aic.sim[aic.sim$actualDist == a & aic.sim$prop==pr & aic.sim$fittedDist == 'tnegb', 'daic']), 
                 col=cols['tnegb'], yaxt='n', xaxt='n', main='', xlim=xlim, lwd=3, 
                 panel.first=rect(par('usr')[1], par('usr')[3], 2, par('usr')[4], col='gray', border=NA))
        }
        
        abline(v=0, col=cols[a], lwd=2)
        
        if(pr == 1) mtext(switch(a, 'fish'='Logseries',
                                 'plnorm'='PoisLogNorm',
                                 'stick'='BrokenStick',
                                 'tnegb'='TruncNegBin'), 
                          side=3, line=1, cex=0.8, col=cols[a])
        if(pr == 4) axis(1, at=aat)
    }
}

mtext(expression(Delta*'AIC'), side=1, line=2, outer=TRUE)
mtext('Relative density', side=2, line=1, outer=TRUE)

dev.off()