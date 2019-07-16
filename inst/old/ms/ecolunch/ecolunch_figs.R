setwd('~/Dropbox/Research/happySAD/ms/ecolunch/figs')
library(socorro)
library(meteR)
library(pika)
library(plyr)

## load hawaii METE objects
load('../../../../hawaiiMETE/mete.byTS.RData')

ylim <- c(1, 10^ceiling(log(max(sapply(mete.byTS, function(x) max(x$sad$data))), 10)))

par(mfcol = c(3, 5), mar = c(0.75, 0.2, 0.75, 0.2), oma = c(2, 2, 0, 0) + 0.5, 
    mgp = c(2, 0.75, 0))

for(i in rev(c('VO', 'LA', 'KH', 'MO', 'KA'))) {
    for(j in c('P', 'H', 'D')) {
        thisSAD <- mete.byTS[[which(byTrophBySite$Site == i & 
                                         byTrophBySite$trophic == j)]]$sad
        plot(thisSAD, 
             ptype = 'rad', log = 'y',
             xaxt = 'n', yaxt = 'n',
             ylim = ylim)
        if(i == 'KA') logAxis(2, expLab = TRUE)
        
        lines(meteDist2Rank(thisSAD), col = hsv(0.6, 0.8, 0.95), lwd = 2)
        
        tnegb <- sad2Rank(sad(thisSAD$data, model = 'plnorm'))
        lines(tnegb, col = hsv(0.1), lwd = 2)
    }
}

