setwd('~/Dropbox/Research/happySAD/ms/esa2016/figs')
devtools::load_all('../../../../pika') 
devtools::load_all('../../../../socorro')
library(maps)
library(mapdata)
library(plyr)

## BBS example
bbsYear <- 2009
bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', bbsYear), as.is = TRUE)
this.sad <- sad(bbs$abund[bbs$route==unique(bbs$route)[2]], 'tnegb', keepData = TRUE)

pdf('ms/esa2016/figs/fig_eg-sad.pdf', width = 3, height = 3)
par(mar=c(2, 2, 0, 0) + 0.1, mgp=c(1, 1, 0))

plot(sort(this.sad$data, TRUE), log = 'y', axes = FALSE, 
     xlab='Species', ylab='Abundance', ylim = c(1, 100))
axis(1, at=c(0, 50), labels = NA)
logAxis(2, labels = NA)

lines(sad2Rank(this.sad), lwd = 2, col = hsv(0, 0.7, 0.9))

dev.off()

## BBS map

bbsInfo <- read.csv('~/Research/datasets/bbs/db/bbsRouteInfo.csv', as.is = TRUE)
bbsSAD <- read.csv('../../../bbsSAD.csv', as.is = TRUE)
bbsSAD <- cbind(bbsSAD, bbsInfo[match(bbsSAD$route, bbsInfo$route.id), c('Longi', 'Lati')])
rownames(bbsSAD) <- NULL


## function to return another function that maps a quantitative value to a color
## after transformation by the function xfun
colVal <- function(x, xfun, cols = hsv(c(0.55, 0.75), c(0.3, 0.8), c(1, 0.5))) {
    cfun <- colorRamp(cols)
    
    newx <- xfun(x)
    m <- 1 / diff(range(newx))
    
    function(v) rgb(cfun(m * (xfun(v) - min(newx))), maxColorValue = 255)
}


## negb z values 

pdf('fig_bbsNegbZ.pdf', width=5, height=5)

with(bbsSAD[!is.na(bbsSAD$z), ], {
    zfun <- colVal(z, function(x) -log(x, 10))
    zrange <- range(-log(z, 10))
    
    d <- density(-log(z, 10), from = zrange[1], to = zrange[2])
    
    layout(matrix(1:2, ncol = 1), heights = c(1, 3))
    par(mar=c(1, 2, 0, 2) + 0.1, mgp = c(1.75, 0.75, 0))
    
    plot(d, axes = FALSE, 
         xlim = rev(zrange),
         ylim = c(min(d$y) - 0.3*diff(range(d$y)), max(d$y)), 
         main = '', xlab = '', ylab = '')
    mtext(expression(z^2*'-value'), side=1, line=2)
    
    rect(xleft = seq(zrange[1], zrange[2] - 0.01*diff(zrange), by = 0.01*diff(zrange)), 
         xright = seq(zrange[1] + 0.01*diff(zrange), zrange[2], by = 0.01*diff(zrange)),
         ybottom = par('usr')[3], ytop = min(d$y) - 0.01*diff(range(d$y)), 
         col = zfun(10^-seq(zrange[1], zrange[2] - 0.01*diff(zrange), by = 0.01*diff(zrange))),
         border = NA)
    
    axis(1, at = seq(0, 10, by = 2), labels = 10^-seq(0, 10, by = 2))
    
    par(mar=c(0.1, 2.1, 0, 2.1))
    
    map('world', c('usa', 'canada'), xlim = c(-180, -50), col=NA)
    map('world', c('canada', 'alaska'), add=TRUE, col='gray', fill=TRUE, border='gray')
    map('world', 'USA:alaska', add=TRUE, col='gray', fill=TRUE, border='gray')
    map('usa', add=TRUE, col='gray', fill=TRUE, border='gray')

    points(Longi, Lati, col=zfun(z), pch=16, cex=0.4)
})

dev.off()


## map of AIC winner
bbsAICWin <- ddply(bbsSAD, 'route', function(x) {
    if(x$dAIC[x$mod == 'tnegb'] <= 2) m <- 'tnegb'
    else m <- x$mod[x$dAIC == 0]
    
    return(data.frame(Longi=x$Longi[1], Lati=x$Lati[1], mod=m))
})

pdf('fig_bbsAIC.pdf', width=5, height=5)

modCol <- hsv(c(0.52, 0.02, 0.12, 0.65), c(0.8, 0.7, 0.8, 0.7), c(0.8, 0.75, 0.9, 0.6))
names(modCol) <- c('fish', 'plnorm', 'stick', 'tnegb')

with(bbsAICWin, {
    win <- sort(table(mod))
    
    layout(matrix(1:2, ncol = 1), heights = c(1, 3))
    par(mar=c(0, 2, 4, 2) + 0.1, mgp = c(1.75, 0.75, 0))
    
    plot(1, xlim = c(0, sum(win)), type = 'n', xaxs = 'i', axes = FALSE, xlab = '', ylab = '')
    
    rect(xleft = c(0, cumsum(win)[-length(win)]), xright = c(cumsum(win)),
         ybottom = par('usr')[3], ytop = par('usr')[4], 
         col = modCol[names(win)],
         border = NA)
    text(cumsum(win) - win/2, mean(par('usr')[3:4]), labels = win, col = 'white', cex = 1.2)
    
    box()
    
    mtext(c('Pois\nlognorm', 'log-\nseries', 'neg binom'), line = 0.5, at = cumsum(win) - win/2,
          col = modCol[names(win)])
    
    par(mar=c(0.1, 2.1, 0, 2.1))
    
    map('world', c('usa', 'canada'), xlim = c(-180, -50), col=NA)
    map('world', c('canada', 'alaska'), add=TRUE, col='gray', fill=TRUE, border='gray')
    map('world', 'USA:alaska', add=TRUE, col='gray', fill=TRUE, border='gray')
    map('usa', add=TRUE, col='gray', fill=TRUE, border='gray')

    points(Longi, Lati, col=modCol[as.character(mod)], pch=16, cex=0.4)
})

dev.off()


## negb k values 

bbsNegB <- ddply(bbs, 'route', function(x) sad(x$abund, model = 'tnegb')$MLE)
bbsNegB <- cbind(bbsNegB, bbsInfo[match(bbsNegB$route, bbsInfo$route.id), c('Longi', 'Lati')])
names(bbsNegB) <- c('route', 'mu', 'k', 'Longi', 'Lati')
rownames(bbsNegB) <- NULL


pdf('fig_bbsNegbK.pdf', width=5, height=5)

with(bbsNegB, {
    zfun <- colVal(k, function(x) x, cols = hsv(c(0.15, 0.03), c(0.55, 0.8), c(0.95, 0.5)))
    zrange <- range(k)
    
    d <- density(k, from = zrange[1], to = zrange[2])
    
    layout(matrix(1:2, ncol = 1), heights = c(1, 3))
    par(mar=c(1, 2, 0, 2) + 0.1, mgp = c(1.75, 0.75, 0))
    
    plot(1, axes = FALSE, type = 'n',
         xlim = zrange, 
         ylim = c(min(d$y) - 0.3*diff(range(d$y)), max(d$y)), 
         main = '', xlab = '', ylab = '')
    mtext('Neg binom k parameter', side=1, line=2)
    axis(1)
    
    rect(xleft = seq(zrange[1], zrange[2] - 0.01*diff(zrange), by = 0.01*diff(zrange)), 
         xright = seq(zrange[1] + 0.01*diff(zrange), zrange[2], by = 0.01*diff(zrange)),
         ybottom = par('usr')[3], ytop = min(d$y) - 0.01*diff(range(d$y)), 
         col = zfun(seq(zrange[1], zrange[2] - 0.01*diff(zrange), by = 0.01*diff(zrange))),
         border = NA)
    rect(xleft = zrange[1], xright = zrange[2],
         ybottom = par('usr')[3], ytop = min(d$y) - 0.01*diff(range(d$y)))
    
    par(mar=c(0.1, 2.1, 0, 2.1))
    
    map('world', c('usa', 'canada'), xlim = c(-180, -50), col=NA)
    map('world', c('canada', 'alaska'), add=TRUE, col='gray', fill=TRUE, border='gray')
    map('world', 'USA:alaska', add=TRUE, col='gray', fill=TRUE, border='gray')
    map('usa', add=TRUE, col='gray', fill=TRUE, border='gray')

    points(Longi, Lati, col=zfun(k), pch=16, cex=0.4)
})

dev.off()


## lognorm and gamma examples
pdf('fig_lnorm.pdf', width = 3.25, height = 3.25)
par(mar=c(3, 3, 0.2, 0.2) + 0.1, mgp=c(2, 0.75, 0), cex.lab=1.2)
curve(dlnorm(x, meanlog = 1.5, sdlog = 0.4), from = 0, to = 14, 
      lwd=4, col='gray45', ylim=c(0, 0.25),
      xlab = 'Latent n', ylab='Probability density')
dev.off()

pdf('fig_gamma.pdf', width = 3.25, height = 3.25)
par(mar=c(3, 3, 0.2, 0.2) + 0.1, mgp=c(2, 0.75, 0), cex.lab=1.2)
curve(dgamma(x, shape = 1.5, scale = 2), from = 0, to = 14, 
      lwd=4, col='gray45', ylim=c(0, 0.25),
      xlab = 'Latent n', ylab='Probability density')
dev.off()



## negbinom example

pdf('fig_negb.pdf', width = 6.5/1.5, height = 5/1.5)
par(mar=c(3, 3, 0, 5) + 0.1, mgp=c(2, 0.75, 0), cex.lab=1.2)

plot(dtnegb(1:14, 8, 0.01), type='b', col=hsv(0.15, 1, 1), lwd=2, cex=1.5, 
     ylim = c(dtnegb(1, 9.2, 100), dtnegb(1, 8, 0.01)), 
     xlab='n', ylab='Probability', axes = FALSE, frame.plot = TRUE)
axis(1)
axis(2)
#, at = seq(0, 0.12, by = 0.04))
points(dtnegb(1:14, 8, 1), type='b', col=hsv(0.35, 0.9, 0.7), lwd=2, cex=1.5)
points(dtnegb(1:14, 9.2, 100), type='b', col=hsv(0.55, 0.5, 0.5), lwd=2, cex=1.5)

par(xpd = NA)
text(14.75, dtnegb(14, 8, 0.01), labels = 'k = 0.01', adj=0, cex=1.5, col=hsv(0.15, 1, 1))
text(14.75, dtnegb(14, 8, 1), labels = 'k = 1', adj=0, cex=1.5, col=hsv(0.35, 0.9, 0.7))
text(14.75, dtnegb(14, 9.2, 100), labels = 'k = 100', adj=0, cex=1.5, col=hsv(0.55, 0.5, 0.5))
par(xpd = FALSE)

dev.off()

