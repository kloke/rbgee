\name{diaggee}
\alias{diaggee}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{diaggee
%%  ~~function to do ... ~~
}
\description{Computes the diagnostics TDBETAS and CFITS between the Wilcoxon and HBR GEE fits.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diaggee(betaw, betahbr, vcw, hess)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betaw}{Coefficients of Wilcoxon fit.
%%     ~~Describe \code{betaw} here~~
}
  \item{betahbr}{Coefficients of HBR fit.
%%     ~~Describe \code{betahbr} here~~
}
  \item{vcw}{Wilcoxon variance-covariance matrix.
%%     ~~Describe \code{vcw} here~~
}
  \item{hess}{Hession.
%%     ~~Describe \code{hess} here~~
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
function (betaw, betahbr, vcw, hess) 
{
    p <- length(betaw)
    n <- length(hess[, 1])
    diff <- betaw - betahbr
    tdbeta = t(cbind(diff)) \%*\% solve(vcw) \%*\% cbind(diff)
    bmtd = (4 * (p + 1)^2)/n
    diffc = hess \%*\% diff
    diffvc = hess \%*\% vcw \%*\% t(hess)
    cfit = diffc/(sqrt(diag(diffvc)))
    bmcf = 2 * sqrt((p + 1)/n)
    list(tdbeta = c(tdbeta), bmtd = bmtd, cfit = c(cfit), bmcf = bmcf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
