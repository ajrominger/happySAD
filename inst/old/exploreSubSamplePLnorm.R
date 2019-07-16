library(socorro)
library(pika)

S <- 1000
n <- 2000
x <- rlnorm(S, 15, 3)
d <- density(log(x))
d$x <- exp(d$x)
plot(d, log = 'x', xaxt = 'n')
logAxis(1, expLab = TRUE)

comm <- sample(1:S, size = n, replace = TRUE, prob = x)
comm <- sort(as.numeric(table(comm)), TRUE)
plot(comm, log = 'y', yaxt = 'n')
logAxis(2)

thisSAD <- sad(comm, model = 'plnorm', keepData = TRUE)
plot(thisSAD, ptype = 'rad', log = 'y')

fitSAD(comm, c('fish', 'plnorm'), keepData = TRUE)
