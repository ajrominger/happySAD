setwd('~/Dropbox/Research/happySAD')
devtools::load_all('../pika')


## important simulation info
# nsim <- 500
# nrep <- 500
# nspp <- 100
# tpois.la <- 4

load('sim_out_Z.RData')

simZ <- array(unlist(sim.outZ), dim=c(dim(sim.outZ[[1]]), length(sim.outZ)),
              dimnames=list(rownames(sim.outZ[[1]]), colnames(sim.outZ[[1]]), 1:length(sim.outZ)))

## flip pval for ll
# simZ[3, 1, ] <- 1 - simZ[3, 1, ] + 1/500

## plot AIC
plot(density(simZ[1, 1, ]), xlim=range(simZ[1, , ]))
lines(density(simZ[1, 2, ]), col='red')


## ==============================
## plot densities of p and z vals
## ==============================

par(mfcol=c(7, 2), mar=rep(0.1, 4), oma=c(4, 4, 1, 1)+0.1)

## plot p-vals
for(i in seq(3, 15, by=2)) {
    plot(density(simZ[i, 1, ]), xlim=range(simZ[seq(3, 15, by=2), , ]), axes=FALSE, frame.plot=TRUE, type='l', main='')
    lines(density(simZ[i, 2, ]), col='red')
    legend('topright', legend=rownames(simZ)[i])
}

## plot z-vals
for(i in seq(2, 14, by=2)) {
    plot(density(log(simZ[i, 1, ])), xlim=log(range(simZ[seq(2, 14, by=2), , ])), axes=FALSE, frame.plot=TRUE, main='')
    lines(density(log(simZ[i, 2, ])), col='red')
}


## ============================================
## plot correlation of p and z vals to deltaAIC
## ============================================

par(mfcol=c(7, 2), mar=rep(0.1, 4), oma=c(4, 4, 1, 1)+0.1)

## plot p-vals
for(i in seq(3, 15, by=2)) {
    plot(simZ[1, 1, ] - simZ[1, 2, ], simZ[i, 1, ], ylim=range(simZ[seq(3, 15, by=2), 1, ]))
}

## plot z-vals
for(i in seq(2, 14, by=2)) {
    plot(simZ[1, 1, ] - simZ[1, 2, ], simZ[i, 1, ], ylim=range(simZ[seq(2, 14, by=2), 1, ]))
}


for(i in seq(2, 14, by=2)) {
    plot(simZ[1, 1, ] - simZ[1, 2, ], simZ[i, 1, ])
}


for(i in seq(2, 14, by=2)) {
    plot(simZ[1, 2, ], simZ[i, 1, ], ylim=range(simZ[seq(2, 14, by=2), 1, ]))
}

for(i in seq(2, 14, by=2)) {
    plot(simZ[1, 2, ], simZ[i, 1, ])
}