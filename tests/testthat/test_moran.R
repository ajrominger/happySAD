# library(ape)
# 
# 
# tr <- rtree(30)
# x <- rnorm(30)
# w <- 1/cophenetic(tr)
# diag(w) <- 0
# Moran.I(x, w, alternative = 'two.sided')
# 
# moranInvDist(x, w, alternative = 'two.sided', B = 999)
