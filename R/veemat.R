veemat <-
function (ehat, center, structure = "CS", substructure = "MM", 
    scores) 
{
    numc <- max(center)
    vc <- c()
    n <- length(ehat)
    v12 <- matrix(rep(0, n^2), ncol = n)
    v12inv <- matrix(rep(0, n^2), ncol = n)
    if (structure == "CS") {
        if (substructure == "DHL") {
            vc <- veedhl(ehat, center)$vc
        }
        if (substructure == "MM") {
            vc <- veemm(ehat, center)$vc
        }
        i1 <- 1
        i2 <- 1
        for (i in 1:numc) {
            esee = ehat[center == i]
            ni <- length(esee)
            i2 <- i1 + ni - 1
            vt <- vc[1] * ((1 - vc[4]) * diag(rep(1, ni)) + vc[4] * 
                rep(1, ni) %*% t(rep(1, ni)))
            vt12 <- sigma12(vt)
            vt12inv <- sigma12inv(vt)
            v12[i1:i2, i1:i2] <- vt12
            v12inv[i1:i2, i1:i2] <- vt12inv
            i1 <- i1 + ni
        }
    }
    if (structure == "WI") {
        v12 <- diag(rep(1, n))
        v12inv <- diag(rep(1, n))
    }
    if (structure == "AR") {
        vc <- veear1m(ehat, center, scores)
        i1 <- 1
        i2 <- 1
        for (i in 1:numc) {
            esee = ehat[center == i]
            ni <- length(esee)
            i2 <- i1 + ni - 1
            kap <- vc[1]/(1 - vc[2]^2)
            vt <- matrix(rep(0, ni^2), ncol = ni)
            for (ip in 1:ni) {
                for (jp in 1:ni) {
                  vt[ip, jp] <- kap * vc[2]^abs((jp - ip))
                }
            }
            vt12 <- sigma12(vt)
            vt12inv <- sigma12inv(vt)
            v12[i1:i2, i1:i2] <- vt12
            v12inv[i1:i2, i1:i2] <- vt12inv
            i1 <- i1 + ni
        }
    }
    list(v12 = v12, v12inv = v12inv, vc = vc)
}
