veemm <-
function (ehat, center) 
{
    numc <- max(center)
    meds <- c(0)
    ehatp <- ehat
    n <- length(ehat)
    for (i in 1:numc) {
        ehatc <- ehat[center == i]
        meds[i] <- median(ehatc)
    }
    reff <- meds
    sigb2 <- mad(meds)^2
    xmat <- model.matrix(~as.factor(center) - 1)
    ehatp <- ehat - xmat %*% meds
    sige2 <- mad(ehatp)^2
    sigt2 <- sigb2 + sige2
    rho <- sigb2/sigt2
    vc <- c(sigt2, sigb2, sige2, rho)
    list(vc = vc, reff = reff)
}
