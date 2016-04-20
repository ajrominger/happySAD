setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(gambin)

load('gambin_sim.RData')

gb.sim <- cbind(model=rownames(gb.sim), as.data.frame(gb.sim))

llComp <- with(gb.sim[gb.sim$model != 'fish', ], aggregate(list(ll=ll), 
                                                           list(truPar=truPar, model=model), 
                                                           function(x) c(mean=mean(x), sd=sd(x))))

plot(llComp$truPar, llComp$ll[, 'mean'], log='x')
