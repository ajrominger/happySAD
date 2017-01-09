setwd('~/Dropbox/Research/happySAD/ms/esa2016/figs')
devtools::load_all('../../../../pika') 
devtools::load_all('../../../../socorro')
library(plyr)
library(parallel)


## BBS 
bbsInfo <- read.csv('~/Research/datasets/bbs/db/bbsRouteInfo.csv', as.is = TRUE)
birdInfo <- read.csv('~/Research/datasets/bbs/db/speciesTableBody.csv', as.is = TRUE)
bbsYear <- 1967:2011

## make site by species by year array
bbsSiteSpp <- mclapply(bbsYear, mc.cores = 6, FUN = function(i) {
    bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', i), as.is = TRUE)
    bbs$route <- factor(bbs$route, levels = sort(unique(bbsInfo$route.id)))
    bbs$spp <- factor(bbs$spp, levels = sort(unique(birdInfo$sppKey)))
    
    ## make site by spp matrix
    out <- tidy2mat(bbs$route, bbs$spp, bbs$abund)
    
    ## make routes not run that year be NA
    out[rowSums(out) == 0, ] <- NA
    
    return(out)
})

bbsSiteSpp <- array(unlist(bbsSiteSpp), dim = c(dim(bbsSiteSpp[[1]]), length(bbsYear)))

## get means across years for each spp at each site (mean is only needed param for Poisson,
## thus we're estimating it under assumption that abundances are approx pois across years)
bbsSiteSppMean <- mclapply(1:dim(bbsSiteSpp)[1], mc.cores = 6, 
                           FUN = function(i) {
                               out <- rowMeans(bbsSiteSpp[i, , ], na.rm = TRUE)
                               out[out == 0] <- NA
                               return(out)
                           })

bbsSiteSppMean <- do.call(rbind, bbsSiteSppMean)


## test poisson-gamma mixture leading to nbinom

i <- 1

## fit gamma
param <- MASS::fitdistr(bbsSiteSppMean[i, ][!is.na(bbsSiteSppMean[i, ])], 'gamma')

sppYear <- bbsSiteSpp[i, , ]
sppYear <- sppYear[rowSums(sppYear, na.rm = TRUE) > 0, !is.na(sppYear[i, ])]

plot(simpECDF(sppYear[, 1]), type = 'l', ylim = c(0.4, 1))
for(j in 2:ncol(sppYear)) lines(simpECDF(sppYear[, j]))

## plot neg binom with param from gamma
curve(pnbinom(x, param$estimate[1], param$estimate[2] / (1 + param$estimate[2])), 
      col = 'red', add = TRUE)

i <- i + 1

