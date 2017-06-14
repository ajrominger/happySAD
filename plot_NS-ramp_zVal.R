library(pika)

setwd('~/Dropbox/Research/happySAD')

allZ <- read.csv('sim_NS-ramp_zVal.csv', as.is = TRUE)

plot(allZ$s, allZ$k1)
abline(h = qchisq(0.95, 1), col = 'red')
