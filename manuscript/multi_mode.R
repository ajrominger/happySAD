library(pika)

preston <- function(x) {

}

internalModes <- function(x, minx = 1) {
    x <- c(x[1], x, minx)
}

x <- c(4, 1, 1, 2, 3, 1)
plot(0:length(x), c(0, x), ylim = c(-max(x), max(x)))
arrows(x0 = 1:length(x), y0 = 0, y1 = diff(c(0, x)), col = 'red')
y <- diff(c(0, x))
rle(y > 0)$lengths

localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    # rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if(length(x) > 2 & x[1] == x[2] & x[2] <= x[3]) {
        y <- y[-1]
    }

    y
}

reduceLocalMaxima <- function(x, by = 1) {
    ii <- localMaxima(x)

}

.expandii <- function(x, ii) {
    newii <- outer(ii, -1:1, '+')
    lapply(ii, function(i) {
        newi <- i + (-length(x):length(x))

    })
}

.expandii(1, 1:4)

x <- c(4, 4, 1, 2, 3, 1)
ii <- localMaxima(x)
x[ii] <- x[ii] - 1
localMaxima(x)


find_peaks <- function (x, by = 1, xmin = 0){
    m <- 1
    x <- c(0, x, 0)

    shape <- diff(sign(diff(x, na.pad = FALSE)))

    pks <- lapply(which(shape < 0), function(i) {
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)

        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))

        if(any(x[c(z:i, (i + 2):w)] <= x[i + 1] - by)) {
            return(i + 1)
        } else {
            return(numeric(0))
        }
    })

    pks <- unlist(pks)

    return(pks - 1)
}


x <- c(3, 4, 1, 2, 2, 1, 1, 0)
plot(rle(x)$values)
find_peaks(rle(x)$values, by = 1)

plot(x)
findpeaks(x, npeaks = 1)





x <- c(1,2,9,9,2,1,1,5,5,1,3)
ii <- localMaxima(x) # 3, 8
plot(1:length(x), x)
points((1:length(x))[ii], x[ii], col = 'red')

x <- c(2,2,9,9,2,1,1,5,5,1)
localMaxima(x) # 3, 8
x <- c(3,2,9,9,2,1,1,5,5,1)
localMaxima(x) # 1, 3, 8



bb <- seq(0.005, 0.1, length.out = 3)
n <- 500
nrep <- 1000

ii <- rep(1:nrep, n)
x <- rfish(n * nrep, bb[1])
lapply(1:nrep, function(i) {
    thisX <- x[ii == i]

})

x <- c(4:1, 4)
plot(x)
u <- l <- matrix(0, nrow = length(x), ncol = length(x))
u[upper.tri(u)] <- 1
l[upper.tri(l)] <- 1

u %*% sign(diff(c(x[1], x)))
l %*% rev(diff(rev(c(x[1],x))))

x <- c(x[1], x, x[length(x)])


plot(x)
x <- c(x[1], x, x[length(x) - 1])

plot(x, col = c('black', c('black', 'red')[as.integer(diff(diff(x) >= 0) < 0) + 1]))


x <- cumsum(cumsum(1:10))
diff(x, lag = 2)
diff(x, differences = 2)
