setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')


## important simulation info
# nsim <- 500
# nrep <- 500
# nspp <- 100
# tpois.la <- 4

load('sim_out_Z.RData')

bla <- sim.outZ[1:3]

array(unlist(bla), dim=c(dim(bla[[1]]), length(bla)))

array(sim.outZ, dim=c(dim(sim.outZ[[1]]), length(sim.outZ)))[1, , ]
