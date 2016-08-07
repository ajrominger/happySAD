setwd('~/Dropbox/Research/pika/R')
sapply(list.files(), source)
devtools::load_all('~/Dropbox/Research/socorro')

setwd('~/Dropbox/Research/happySAD/ms/esa2016/figs')

pdf('fig_negb.pdf', width = 6.5/1.5, height = 5/1.5)
par(mar=c(3, 3, 0, 5) + 0.1, mgp=c(2, 0.75, 0), cex.lab=1.2)

plot(dtnegb(1:14, 8, 0.01), type='b', col=hsv(0.15, 1, 1), lwd=2, cex=1.5, 
     ylim = c(dtnegb(1, 9.2, 100), dtnegb(1, 8, 0.01)), 
     xlab='n', ylab='Probability', axes = FALSE, frame.plot = TRUE)
axis(1)
axis(2)
#, at = seq(0, 0.12, by = 0.04))
points(dtnegb(1:14, 8, 1), type='b', col=hsv(0.35, 0.9, 0.7), lwd=2, cex=1.5)
points(dtnegb(1:14, 9.2, 100), type='b', col=hsv(0.55, 0.5, 0.5), lwd=2, cex=1.5)

par(xpd = NA)
text(14.75, dtnegb(14, 8, 0.01), labels = 'k = 0.01', adj=0, cex=1.5, col=hsv(0.15, 1, 1))
text(14.75, dtnegb(14, 8, 1), labels = 'k = 1', adj=0, cex=1.5, col=hsv(0.35, 0.9, 0.7))
text(14.75, dtnegb(14, 9.2, 100), labels = 'k = 100', adj=0, cex=1.5, col=hsv(0.55, 0.5, 0.5))
par(xpd = FALSE)

dev.off()
