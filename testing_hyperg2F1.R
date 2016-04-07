## ============================================
## testing bounds of convergance for hyperg_2F1
## ============================================

p <- seq(0.001, 0.0008, length=80)
z <- exp(-p)
a <- floor(seq(1, 20000, length=300))

y <- sapply(a, function(A) gsl::hyperg_2F1(A, 1, A+1, z))
y[is.nan(y)] <- NA
image(z, a, y)

bound <- apply(y, 1, function(x) {
    thisone <- which(is.na(x))[1]
    if(length(thisone) > 0) return(thisone)
    else return(NA)
})

zcrit <- z[!is.na(bound)]
acrit <- a[bound][!is.na(bound)]

points(zcrit, acrit)
curve(exp(-34)*x^-46653, add=TRUE)

