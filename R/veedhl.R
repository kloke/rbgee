veedhl <-
function (ehat, center) 
{
    numc <- max(center)
    hls <- c(0)
    ehatp <- ehat
    n <- length(ehat)
    for (i in 1:numc) {
        ehatc <- ehat[center == i]
        hls[i] <- hlest(ehatc)
    }
    reff <- hls
    sigb2 <- (pi/3) * ((disp(0, hls, hls, wscores)/numc)^2)
    xmat <- model.matrix(~as.factor(center) - 1)
    ehatp <- ehat - xmat %*% hls
    sige2 <- (pi/3) * ((disp(0, ehatp, ehatp, wscores)/n)^2)
    sigt2 <- sigb2 + sige2
    rho <- sigb2/sigt2
    vc <- c(sigt2, sigb2, sige2, rho)
    list(vc = vc, reff = reff)
}
