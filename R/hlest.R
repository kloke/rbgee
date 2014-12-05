hlest <-
function (x) 
{
    mat <- pairup(x, type = "leq")
    hlest <- median((mat[, 1] + mat[, 2])/2)
    return(hlest)
}
