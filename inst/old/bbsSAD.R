setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika') 
devtools::load_all('../socorro')
library(plyr)


## read BBS

bbsYear <- 2009
bbs <- read.csv(sprintf('~/Research/datasets/bbs/db/bbs%s.csv', bbsYear), as.is = TRUE)


## fit all SADs to BBS data

bbsSAD <- ddply(bbs, 'route', function(x) {
    fit <- fitSAD(x$abund, models = c('tnegb', 'plnorm', 'fish'), keepData = TRUE)
    aic <- sapply(fit, AIC)
    z <- rep(NA, length(aic))
    z[1] <- logLikZ(fit$tnegb)$z
    
    return(data.frame(mod=names(fit), aic=aic, dAIC=aic-min(aic), z=z))
})

## write output
write.table(bbsSAD, file = 'bbsSAD.csv', sep = ',', row.names = FALSE)
