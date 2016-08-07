## loading pika without dependencies
# devtools::load_all('../pika') 
setwd('~/Dropbox/Research/pika/R')
sapply(list.files(), source)
devtools::load_all('~/Dropbox/Research/socorro')
library(plyr)

setwd('~/Dropbox/Research/happySAD')

## read BBS
bbsYear <- 2009
bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', bbsYear), as.is = TRUE)

## read BBS info
bbsInfo <- read.csv('~/Research/datasets/bbs/db/bbsRouteInfo.csv', as.is = TRUE)


## make a quick plot for ESA presentation
this.sad <- sad(bbs$abund[bbs$route==unique(bbs$route)[2]], 'tnegb', keepData = TRUE)

pdf('ms/esa2016/figs/fig_eg-sad.pdf', width = 3, height = 3)
par(mar=c(2, 2, 0, 0) + 0.1, mgp=c(1, 1, 0))

plot(sort(this.sad$data, TRUE), log = 'y', axes = FALSE, 
     xlab='Species', ylab='Abundance', ylim = c(1, 100))
axis(1, at=c(0, 50), labels = NA)
logAxis(2, labels = NA)

lines(sad2Rank(this.sad), lwd = 2, col = hsv(0, 0.7, 0.9))

dev.off()


## ============================================
## the real stuff: fitting all SADs to BBS data
## ============================================

bbsSAD <- ddply(bbs[bbs$route %in% unique(bbs$route)[1:10], ], 'route', function(x) {
    fit <- fitSAD(x$abund, models = c('tnegb', 'fish'), keepData = TRUE)
    aic <- sapply(fit, AIC)
    z <- rep(NA, length(aic))
    z[1] <- logLikZ(fit$tnegb)$z
    
    return(data.frame(mod=names(fit), aic=aic, dAIC=aic-min(aic), z=z))
})

