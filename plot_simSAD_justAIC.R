## set up
setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')
library(reshape2)

## read in files
files <- list.files('.', pattern='sim_out_aic')
sim.aic <- c()
for(f in files) {
    load(f)
    sim.aic <- c(sim.aic, sim.out.aic)
}

## clean simulation
aic.sim <- melt(sim.aic)
colnames(aic.sim) <- c('stat', 'fittedDist', 'value', 'pars', 'actualDist', 'prop', 'nspp', 'rep')
aic.sim <- dcast(aic.sim, actualDist + pars + fittedDist + nspp + prop + rep ~ stat, value.var='value')

