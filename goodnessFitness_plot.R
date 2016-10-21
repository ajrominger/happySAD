library(pika)
library(socorro)
library(plyr)

setwd('~/Dropbox/Research/happySAD')

goodness <- read.csv('sim_goodnessFitness.csv', as.is = TRUE)