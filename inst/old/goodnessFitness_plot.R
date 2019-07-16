library(pika)
library(socorro)
library(plyr)
library(reshape2)

setwd('~/Dropbox/Research/happySAD')

## read in simulation results
goodness <- read.csv('goodnessFitness_res.csv', as.is = TRUE)
goodness <- melt(goodness, id.vars = c('fg', 'ff', 'S'))

## summarize
goodSumm <- ddply(goodness, c('fg', 'ff', 'S', 'variable'), function(d) {
    out <- c(mean(d$value), quantile(d$value, probs = c(0.025, 0.975)))
    names(out) <- c('mean', 'ciLo', 'ciHi')
    return(out)
})


## plotting

v <- 'z_ll'
maxN <- 30

layout(matrix(1:(2 * (length(fg) * (length(ff) - 1))), ncol = 2, byrow = TRUE))
par(mar = c(0, 2, 0, 1) + 0.5, oma = c(3, 0, 0,  0))

for(i in unique(goodSumm$fg)) {
    for(j in unique(goodSumm$ff)) {
        ## skip plotting when generating and fitted model are the same
        if(i != j) {
            ## PMF for generating function
            g <- dtnegb(1:maxN, as.numeric(gsub(',.*', '', i)), as.numeric(gsub(',.*', '', i)))
            
            ## PMF for fitted function
            f <- dtnegb(1:maxN, as.numeric(gsub(',.*', '', j)), as.numeric(gsub(',.*', '', j)))
            
            ## plot PMFs
            plot(g, ylim = range(f, g),  type = 'l')
            points(f, type = 'l', col = 'red')
            legend('topright', legend = c(i, j), col = c('black', 'red'), lty = 1)
            
            ## plot statistic
            with(goodSumm[goodSumm$variable == v & goodSumm$fg == i & 
                              (goodSumm$ff == j | goodSumm$ff == i), ], 
                 plot(S, mean, pch = 21, bg = 'white', 
                      col = c('black', 'red')[(ff == j) + 1], 
                      ylim = range(goodSumm[goodSumm$variable == v, c('mean', 'ciLo', 'ciHi')]), 
                      xaxt = 'n', log = 'x',
                      panel.first = arrows(x0 = S, y0 = ciLo, y1 = ciHi, 
                                           code = 3, angle = 90, length = 0.05,
                                           col = c('black', 'red')[(ff == j) + 1])))
        }
    }
}
logAxis(1, expLab = TRUE)
levels(goodSumm$variable)
