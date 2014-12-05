geerfit <-
function (y, xmat, center, scores = wscores, geemod = "LM", structure = "CS", 
    substructure = "MM", med = TRUE, varcovst = "var2", maxstp = 50, 
    eps = 1e-05, hbrs = FALSE, delta = 0.8, hparm = 2) 
{
    int = TRUE
    xmat <- centerx(xmat)
    n <- length(y)
    one <- rep(1, n)
    p <- length(xmat[1, ])
    if (int) {
        if (hbrs == FALSE) {
            fit0 <- rfit(y ~ xmat, scores = scores)
        }
        if (hbrs == TRUE) {
            fit0 <- hbrfit(y ~ xmat)
        }
        beta0 <- fit0$coef[2:(p + 1)]
        y <- y - fit0$coef[1]
        alphat <- fit0$coef[1]
    }
    else {
        if (hbrs == FALSE) {
            fit0 <- rfit(y ~ xmat, scores = scores)
            yhatint <- fitted.values(fit0)
        }
        if (hbrs == TRUE) {
            fit0 <- hbrfit(y ~ xmat)
            yhatint <- y - fit0$resid
        }
        beta0 <- solve(t(xmat) %*% xmat) %*% t(xmat) %*% yhatint
    }
    colbeta <- matrix(rep(0, maxstp * 2 * p), ncol = 2 * p)
    istp <- 0
    iflag <- 0
    ic <- 1
    while (ic > 0) {
        apbeta0 <- getAp(xmat, center, beta0, geemod)
        S <- y - apbeta0
        vee <- veemat(S, center, structure, substructure, scores)
        vee12 <- as.matrix(vee$v12)
        vee12inv <- as.matrix(vee$v12inv)
        Dmat <- getD(xmat, center, beta0, geemod)
        eitb <- vee12inv %*% (y - apbeta0)
        Dmatb <- vee12inv %*% Dmat
        wtstuff <- wtmat(Dmatb, eitb, med, scores, hbrs)
        wts <- wtstuff$wts
        wts <- as.vector(wts)
        msbeta0 <- wtstuff$mb * vee12 %*% one
        S <- y - apbeta0 - msbeta0
        veeinv <- vee12inv %*% diag(wts) %*% vee12inv
        part1 <- t(Dmat) %*% veeinv %*% Dmat
        inc <- solve(part1) %*% t(Dmat) %*% veeinv %*% S
        beta1 <- beta0 + inc
        istp <- istp + 1
        colbeta[istp, ] <- c(beta0, beta1)
        chk <- sum(inc^2)/(sum(beta0^2) + 1e-05)
        if (chk < eps) {
            iflag <- 1
            ic <- 0
            betag <- beta1
        }
        else {
            if (istp >= maxstp) {
                ic <- 0
                betag <- beta1
            }
            beta0 <- beta1
        }
    }
    apbetag <- getAp(xmat, center, betag, geemod)
    S <- y - apbetag
    vee <- veemat(S, center, structure, substructure, scores)
    vc <- vee$vc
    vee12 <- as.matrix(vee$v12)
    vee12inv <- as.matrix(vee$v12inv)
    Dmat <- getD(xmat, center, betag, geemod)
    eitbg <- vee12inv %*% (y - apbetag)
    Dmatb <- vee12inv %*% Dmat
    wtstuff <- wtmat(Dmatb, eitbg, med, scores, hbrs)
    wts <- wtstuff$wts
    wts <- as.vector(wts)
    Dmat <- getD(xmat, center, betag, geemod)
    if (varcovst == "var2") {
        vcmat <- varcov2(eitbg, center, scores, vee12inv, Dmat, 
            p, delta, hparm)
    }
    if (varcovst == "var1") {
        vcmat <- varcov1(eitbg, center, scores, vee12inv, Dmat, 
            wts)
    }
    if (varcovst == "var3") {
        vcmat <- varcov3(eitbg, scores, vee12inv, Dmat, p, delta, 
            hparm)
    }
    betahist <- colbeta[1:istp, ]
    se <- sqrt(diag(vcmat))
    tab <- cbind(betag, se, betag/se)
    colnames(tab) <- c("Est", "SE", "t-ratio")
    varcov <- vcmat
    ehat1 <- y - xmat %*% betag
    alphat2 <- median(ehat1) + alphat
    ehat <- ehat1 - median(ehat1)
    yhat <- y - ehat
    list(tab = tab, betahist = betahist, varcov = varcov, stanhess = Dmatb, 
        istp = istp, vc = vc, beta0 = alphat2, residuals = ehat, 
        fitted.values = yhat)
}
