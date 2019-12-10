setwd('~/Dropbox/Research/happySAD/inst/old')
library(pika)
library(gambin)

load('gambin_sim.RData')

gb.sim <- cbind(model=rownames(gb.sim), as.data.frame(gb.sim))

llComp <- with(gb.sim, aggregate(list(ll=ll),
                                 list(truPar=truPar, model=model),
                                 function(x) c(mean=mean(x), sd=sd(x))))

plot(llComp$truPar, llComp$ll[, 'mean'], log='x')

with(gb.sim[gb.sim$model=='binFish', ], plot(truPar, MLE, log = 'xy'))
with(gb.sim[gb.sim$model=='fish', ], points(truPar, MLE, col='red'))
abline(0, 1)
