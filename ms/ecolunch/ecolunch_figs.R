setwd('~/Dropbox/Research/happySAD/ms/ecolunch/figs')
library(socorro)
library(meteR)
library(plyr)

## load hawaii METE objects
load('../../../../hawaiiMETE/mete.byTS.RData')

ylim <- c(1, 10^ceiling(log(max(sapply(mete.byTS, function(x) max(x$sad$data))), 10)))

par(mfcol = c(3, 5), mar = c(0.75, 0.2, 0.75, 0.2), oma = c(2, 2, 0, 0) + 0.5, 
    mgp = c(2, 0.75, 0))

for(i in rev(c('VO', 'LA', 'KH', 'MO', 'KA'))) {
    for(j in c('P', 'H', 'D')) {
        plot(mete.byTS[[which(byTrophBySite$Site == i & byTrophBySite$trophic == j)]]$sad, 
             ptype = 'rad', log = 'y', add.legend = FALSE, th.col = 'black',
             xaxt = 'n', yaxt = 'n', xlab = '', ylab = 'n',
             ylim = ylim)
        if(i == 'KA') logAxis(2, expLab = TRUE)
    }
}

