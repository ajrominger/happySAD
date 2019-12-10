moranInvDist <- function(x, invDist, alternative = 'two.sided', B = 999) {
    r <- rowSums(invDist)
    r[r == 0] <- 1
    
    w <- invDist / r
    W <- sum(w)
    xdiff <- x - mean(x)
    n <- length(x)
    
    obs <- n / W * sum(w * outer(xdiff, xdiff)) / sum(xdiff^2)
    
    if(!is.null(B)) {
        nullSim <- replicate(B, {
            newx <- sample(x)
            return(moranInvDist(newx, invDist, B = NULL))
        })
        nullSim <- c(nullSim, obs)
        
        pl <- mean(obs >= nullSim)
        pg <- mean(obs <= nullSim)
        pt <- 2 * min(pl, pg)
        return(list(obs = obs, 
                    p = switch(alternative,
                        'two.sided' = pt,
                        'less' = pl,
                        'greater' = pg
                    )))
    } else {
        return(obs)
    }
}
