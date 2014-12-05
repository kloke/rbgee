hbrwtsnew <-
function (xmat, y, robdis2 = robdist.hbrfit(xmat), percent = 0.95, 
    intest = ltsreg(xmat, y)$coef) 
{
    xmat = as.matrix(xmat)
    y = as.matrix(y)
    n = dim(xmat)[1]
    p = dim(xmat)[2]
    cut = qchisq(percent, p)
    resids = y - intest[1] - xmat %*% as.matrix(intest[2:(p + 
        1)])
    sigma = mad(resids)
    m = psi(cut/robdis2)
    a = resids/(sigma * m)
    c = (median(a) + 3 * mad(a))^2
    h = sqrt(c)/a
    ans = psi(abs(h))
    ans
}
