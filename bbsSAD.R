## loading pika without dependencies
# devtools::load_all('../pika') 
setwd('~/Dropbox/Research/pika/R')
sapply(list.files(), source)
devtools::load_all('~/Dropbox/Research/socorro') 

setwd('~/Dropbox/Research/happySAD')

## read BBS
bbsYear <- 2009
bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', bbsYear), as.is = TRUE)
head(bbs)

## make a quick plot for ESA presentation
this.sad <- sad(bbs$abund[bbs$route==unique(bbs$route)[2]], 'tnegb', keepData = TRUE)


plot(sort(this.sad$data, TRUE), log='y')
lines(sad2Rank(this.sad))


x <- lapply(unique(bbs$route)[1:2], function(r) {
    fit <- fitSAD(bbs$abund[bbs$route == r], models = c('tnegb', 'fish'), keepData = TRUE)
    aic <- sapply(fit, AIC)
    z <- rep(NA, length(aic))
    z[1] <- logLikZ(fit$tnegb)$z
    
    return(cbind(aic, dAIC=aic-min(aic), z))
})

do.call(AIC, x[[1]])
AIC(x[[1]][[1]])
AIC(x[[1]][[2]])


logLikZ(x[[2]][[1]])
