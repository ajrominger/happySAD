#' @title Find local maxima
#'
#' @description Returns a vector of indeces of the local maxima
#'
#' @details Thanks to stackoverflow user stas g for the starting point solution in
#' https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r
#'
#' @param x vector of sequential numbers within which to look for local maxima
#' @param by minimum value by which a maximum must exceed its neighbors
#' @param xmin the hypothetical minimum value \code{x} might achieve
#'
#' @return a vector of indeces of the local maxima
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

localMaxima <- function(x, by = 1, xmin = 0) {
    x <- c(xmin, rle(x)$values, xmin)
    x <- rle(by * round(x / by))$values


    shape <- diff(sign(diff(x, na.pad = FALSE)))

    pks <- lapply(which(shape < 0), function(i) {
        w <- i + 2
        w <- ifelse(w < length(x), w, length(x))

        if(any(x[c(i, w)] <= x[i + 1] - by)) {
            return(i)
        } else {
            return(numeric(0))
        }
    })

    pks <- unlist(pks)

    return(pks)
}
