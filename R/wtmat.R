wtmat <-
function (dmatb, eitb, med = TRUE, scores = wscores, hbrs = FALSE) 
{
    n <- length(eitb)
    u <- (1:n)/(n + 1)
    eps <- 1e-06
    sc <- getScores(scores, u)
    top <- sc[rank(eitb, ties.method = "first")]
    wts <- rep(0, n)
    if (med) {
        mb <- median(eitb)
        ehat <- eitb - mb
        wts <- top/ehat
        wts[is.infinite(wts)] <- 0
        wts[is.nan(wts)] <- 0
        wts[wts == 0] <- max(wts)
    }
    else {
        scl <- sc[sc < 0]
        mp1 <- length(scl) + 1
        eitbs <- sort(eitb)
        mb <- eitbs[mp1]
        ehat <- eitb - mb
        wts <- top/ehat
        wts[is.infinite(wts)] <- 0
        wts[is.nan(wts)] <- 0
        wts[wts == 0] <- max(wts)
    }
    if (hbrs) {
        wts2 <- hbrwtsnew(dmatb, eitb)
        newwts <- wts * wts2
        newwts[is.infinite(wts)] <- 0
        newwts[is.nan(wts)] <- 0
        newwts[wts == 0] <- max(wts)
        wts <- newwts
    }
    list(wts = wts, mb = mb)
}
