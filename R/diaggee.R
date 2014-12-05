diaggee <-
function (betaw, betahbr, vcw, hess) 
{
    p <- length(betaw)
    n <- length(hess[, 1])
    diff <- betaw - betahbr
    tdbeta = t(cbind(diff)) %*% solve(vcw) %*% cbind(diff)
    bmtd = (4 * (p + 1)^2)/n
    diffc = hess %*% diff
    diffvc = hess %*% vcw %*% t(hess)
    cfit = diffc/(sqrt(diag(diffvc)))
    bmcf = 2 * sqrt((p + 1)/n)
    list(tdbeta = c(tdbeta), bmtd = bmtd, cfit = c(cfit), bmcf = bmcf)
}
