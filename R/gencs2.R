gencs2 <-
function (vecni, vc, beta, xmat, itypee = 1, itypeb = 1) 
{
    numc <- length(vecni)
    n <- sum(vecni)
    center <- rep(0, n)
    ind <- rep(0, n)
    errs <- rep(0, n)
    i1 <- 1
    i2 <- 0
    for (j in 1:numc) {
        if (itypeb == 1) {
            bj <- rnorm(1, 0, sqrt(vc[1]))
        }
        i2 <- i2 + vecni[j]
        ni <- vecni[j]
        ni2 <- floor(ni/2)
        center[i1:i2] <- rep(j, ni)
        if (itypee == 1) {
            ei <- rnorm(ni, 0, sqrt(vc[2]))
        }
        errs[i1:i2] <- ei + bj
        i1 <- i1 + ni
    }
    y <- xmat %*% beta + errs
    list(y = y, center = center)
}
