#' @title Plotting SAD shapes
#'  
#' @description plot the shapes of SADs in parameter space
#' 
# @details 
#' 
#' @param pars a list containing the parameter values of interest
#' @param nn vector of integer abundances for which to plot probabilities
#' @param dfun the density function of interest
#' @param main plot title
#' @param ... additional parameters passed to \code{plot}
#' @author Andy Rominger <ajrominger@@gmail.com>
#' 
#' @export

sadShape <- function(pars, nn, dfun, main, ...) {
    npar1 <- length(pars[[1]])
    if(length(pars) == 1) {
        addY <- FALSE
        pars <- c(pars, list(1))
    } else {
        addY <- TRUE
    }
    npar2 <- length(pars[[2]])
    
    
    allpar <- as.matrix(expand.grid(pars))
    
    # browser()
    
    ymax <- max(sapply(1:nrow(allpar), function(i) {
        if(addY) {
            plist <- c(list(nn), as.list(allpar[i, ]))
        } else {
            plist <- c(list(nn), as.list(allpar[i, -2]))
        }
        
        return(max(do.call(dfun, plist)))
    }))
    
    par(plt = c(0.25, 0.9, 0.25, 0.9))
    
    plot(range(allpar[, 1]) + c(-1, 1) * max(diff(allpar[, 1])) / 2, 
         range(allpar[, 2]) + c(-1, 1) * max(diff(allpar[, 2])) / 2,
         type = 'n', yaxt = 'n', axes = FALSE, xaxs = 'i', yaxs = 'i', ...)
    
    axis(1, line = 0.2)
    if(addY) {
        axis(2, line = 0.2)
    }
    
    mtext(main, line = 0.5)
    
    strt1 <- par('plt')[1]
    bystep1 <- diff(par('plt')[1:2]) / npar1
    m1 <- cbind(seq(strt1, by = bystep1, length.out = npar1), 
               seq(strt1 + bystep1, by = bystep1, length.out = npar1))
    
    strt2 <- par('plt')[3]
    bystep2 <- diff(par('plt')[3:4]) / npar2
    m2 <- cbind(seq(strt2, by = bystep2, length.out = npar2), 
                seq(strt2 + bystep2, by = bystep2, length.out = npar2))
    
    foo <- as.matrix(expand.grid(1:npar1, 1:npar2))
    m <- cbind(m1[foo[, 1], ], m2[foo[, 2], ])
    
    split.screen(m, erase = FALSE)
    for(i in 1:nrow(allpar)) {
        screen(i, new = FALSE)
        par(mar = c(0, 0, 0, 0) + 0.1, mgp = c(1, 0.2, 0), tcl = -0.2)
        
        if(addY) {
            plist <- c(list(nn), as.list(allpar[i, ]))
        } else {
            plist <- c(list(nn), as.list(allpar[i, -2]))
        }
        
        y <- do.call(dfun, plist)
        plot(nn, y, xlab = '', ylab = '', log = 'x', 
             type = 'l', lwd = 2, 
             ylim = c(.Machine$double.eps^0.75, ymax), xaxt = 'n', yaxt = 'n')
        polygon(c(nn, nn[1]), c(y, y[length(y)]), col = 'gray')
    }
    
    close.screen(all.screens = TRUE)
}
