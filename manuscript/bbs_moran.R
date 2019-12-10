library(pika)
library(geosphere)
library(ape)
library(sp)
library(gstat)
library(socorro)
library(viridis)

bbs <- read.csv('~/Dropbox/Research/IBL/Data/bbs/bbs2011.csv', as.is = TRUE)
bbsInfo <- read.csv('~/Dropbox/Research/IBL/Data/bbs/bbsRouteInfo.csv', as.is = TRUE)
theseRoutes <- unique(bbs$route)
bbsInfo <- bbsInfo[match(theseRoutes, bbsInfo$route.id), ]

# calculate distance matrix

d <- matrix(0, nrow = length(theseRoutes), ncol = length(theseRoutes))
d[lower.tri(d)] <- 1
routeComp <- which(d == 1, arr.ind = TRUE)

d[lower.tri(d)] <- distGeo(bbsInfo[routeComp[, 1], c('Longi', 'Lati')],
                           bbsInfo[routeComp[, 2], c('Longi', 'Lati')])
d[upper.tri(d)] <- t(d)[upper.tri(d)]


# calculate negbinom stats

bbsByRoute <- split(bbs, bbs$route)

nbInfo <- parallel::mclapply(bbsByRoute, mc.cores = 3, FUN = function(x) {
    s <- sad(x$abund, 'tnegb', keepData = TRUE)
    sfish <- sad(x$abund, 'fish', keepData = TRUE)

    o <- c(s$nobs, sum(x$abund), s$MLE, sfish$MLE,
           logLikZ(s) <= qchisq(0.95, 1), AIC(sfish) - AIC(s))
    names(o) <- c('S', 'N', 'mu', 'k', 'beta', 'nbOK', 'fishDiff')

    return(o)
})

nbInfo <- data.frame(route = names(nbInfo), do.call(rbind, nbInfo))
nbInfo$nbOK <- as.logical(nbInfo$nbOK)
nbInfo$fishBetter <- nbInfo$fishDiff < 0

nbInfo <- nbInfo[match(bbsInfo$route.id, nbInfo$route), ]
nbInfo$lon <- bbsInfo$Longi
nbInfo$lat <- bbsInfo$Lati
coordinates(nbInfo) = ~ lon + lat
proj4string(nbInfo) <- CRS('+proj=longlat +ellps=WGS84')

plot(nbInfo$k, nbInfo$beta, log = 'xy')

plot(nbInfo$mu, nbInfo$k, log = 'xy',
     col = c('black', 'red')[(nbInfo$fishDiff < -2) + 1])

pdf('foo.pdf', width = 14, height = 10)
par(mar = rep(0, 4))
plot(nbInfo, pch = 16, cex = 0.5, col = quantCol(nbInfo$beta, viridis(20), trans = 'log'))
dev.off()

w <- 1 / d
w[!is.finite(w)] <- 0
Moran.I(nbInfo$k, w)

v <- variogram(k ~ 1, data = nbInfo, cloud = TRUE)
v <- do.call(data.frame, unclass(v)[c('dist', 'gamma')])
v <- lapply(split(v, cut(v$dist, 50)), function(x) {
    c(mean(x$dist), mean(x$gamma), quantile(x$gamma, probs = c(0.25, 0.75)))
})
v <- do.call(rbind, v)

nrow(v)
plot(v[, 1:2], ylim = c(0, 0.05))
segments(x0 = v[, 1], y0 = v[, 3], y1 = v[, 4])
