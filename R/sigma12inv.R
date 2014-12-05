sigma12inv <-
function (sigma) 
{
    temp <- eigen(sigma, symmetric = TRUE)
    sigma12 <- temp$vectors %*% diag(temp$values^(-0.5)) %*% 
        t(temp$vectors)
    sigma12
}
