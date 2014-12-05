\name{veear1m}
\alias{veear1m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{veear1m
%%  ~~function to do ... ~~
}
\description{Computes AR(1) estimate for workink covariance matrix.   Structure should be set at "AR".
See Section 8.6.3 of Kloke and McKean (2014).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
veear1m(ehat, center, scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ehat}{
%%     ~~Describe \code{ehat} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
}
  \item{scores}{
%%     ~~Describe \code{scores} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Kloke and McKean (2014), Nonparametrics Using R, Boca Raton: Chapman-Hall.

%% ~put references to the literature/web site here ~
}
\author{Joseph W. McKean
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ehat, center, scores) 
{
    numc <- max(center)
    n <- length(ehat)
    rho <- rep(0, numc)
    sig <- rep(0, numc)
    for (i in 1:numc) {
        ehatc <- ehat[center == i]
        ni <- length(ehatc)
        elp1 <- ehatc[2:ni]
        elm1 <- ehatc[1:(ni - 1)]
        fiti <- rfit(elp1 ~ elm1, scores = scores)
        rho[i] <- fiti$coef[2]
        sig[i] <- mad(fiti$resid)
    }
    rhomed <- median(rho)
    sigmed2 <- median(sig)^2
    if (rhomed > 0.95) {
        rhomed <- 0.95
    }
    if (rhomed < -0.95) {
        rhomed <- -0.95
    }
    vc <- c(sigmed2, rhomed)
    return(vc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
