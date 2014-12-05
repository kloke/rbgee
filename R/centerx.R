centerx <-
function (x) 
{
    x = as.matrix(x)
    n = length(x[, 1])
    one = matrix(rep(1, n), ncol = 1)
    x - (one %*% t(one)/n) %*% x
}
