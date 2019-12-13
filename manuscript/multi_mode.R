library(pika)
library(happySAD)
library(viridis)

s <- seq(100, 500, length.out = 5)
b <- 10^seq(-3, -1, by = 0.5)
sb <- as.matrix(expand.grid(s = s, b = b))

p1mod <- lapply(1:nrow(sb), function(i) {
    thisS <- sb[i, 1]
    thisB <- sb[i, 2]

    dsingModeGivenS(thisS, dfish, thisB)
})


p1modMat <- matrix(unlist(p1mod), nrow = length(s), ncol = length(b))
image(1:length(b), s, t(p1modMat), col = viridis(16), zlim = c(0, 1),
      xlab = expression(beta), ylab = 'S', xaxt = 'n')
axis(1, at = 1:length(b), labels = paste0('10^', log(b, 10)))
