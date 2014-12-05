veear1m <-
function (ehat, center, scores) 
{
    numc <- max(center)
    n <- length(ehat)
    rho <- rep(0, numc)
    sig <- rep(0, numc)
    for (i in 1:numc) {
        ehatc <- ehat[center == i]
        ni <- length(ehatc)
        elp1 <- ehatc[2:ni]
        elm1 <- ehatc[1:(ni - 1)]
        fiti <- rfit(elp1 ~ elm1, scores = scores)
        rho[i] <- fiti$coef[2]
        sig[i] <- mad(fiti$resid)
    }
    rhomed <- median(rho)
    sigmed2 <- median(sig)^2
    if (rhomed > 0.95) {
        rhomed <- 0.95
    }
    if (rhomed < -0.95) {
        rhomed <- -0.95
    }
    vc <- c(sigmed2, rhomed)
    return(vc)
}
