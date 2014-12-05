getAp <-
function (xmat, center, beta, geemod = "LM") 
{
    if (geemod == "LM") {
        apbeta <- xmat %*% beta
    }
    return(apbeta)
}
