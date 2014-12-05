varcov3 <-
function (eitbg, scores, v12inv, Dmat, p, delta, hparm) 
{
    n <- length(eitbg)
    u <- (1:n)/(n + 1)
    sc <- getScores(scores, u)
    tauhat <- gettauF0(eitbg, p, scores, delta, hparm)
    veeinv <- v12inv %*% v12inv
    part1 <- t(Dmat) %*% veeinv %*% Dmat
    varcov3 <- (tauhat^2) * solve(part1)
    return(varcov3)
}
