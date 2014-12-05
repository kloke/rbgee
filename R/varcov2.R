varcov2 <-
function (eitbg, center, scores, v12inv, Dmat, p, delta, hparm) 
{
    numc <- max(center)
    n <- length(eitbg)
    u <- (1:n)/(n + 1)
    sc <- getScores(scores, u)
    rs <- sc[rank(eitbg, ties.method = "first")]
    brrp <- matrix(rep(0, n * n), ncol = n)
    tauhat <- gettauF0(eitbg, p, scores, delta, hparm)
    veeinv <- v12inv %*% v12inv
    part1 <- t(Dmat) %*% veeinv %*% Dmat
    i1 <- 1
    i2 <- 1
    for (i in 1:numc) {
        rsc <- rs[center == i]
        ni <- length(rsc)
        rscbar <- mean(rsc)
        rsc <- rsc - rscbar
        rrp <- rsc %*% t(rsc)
        i2 <- i1 + ni - 1
        brrp[i1:i2, i1:i2] <- rrp
        i1 <- i1 + ni
    }
    mid <- t(Dmat) %*% v12inv %*% brrp %*% v12inv %*% Dmat
    varcov2 <- (tauhat^2) * solve(part1) %*% mid %*% solve(part1)
    return(varcov2)
}
