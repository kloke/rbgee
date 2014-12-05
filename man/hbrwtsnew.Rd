\name{hbrwtsnew}
\alias{hbrwtsnew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{HBR Weights
%%  ~~function to do ... ~~
}
\description{Computes the high breakdown weights used in the HBR GEE estimator
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hbrwtsnew(xmat, y, robdis2 = robdist.hbrfit(xmat), percent = 0.95, intest = ltsreg(xmat, y)$coef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmat}{The design matrix.
%%     ~~Describe \code{xmat} here~~
}
  \item{y}{The response vector.
%%     ~~Describe \code{y} here~~
}
  \item{robdis2}{
%%     ~~Describe \code{robdis2} here~~
}
  \item{percent}{
%%     ~~Describe \code{percent} here~~
}
  \item{intest}{
%%     ~~Describe \code{intest} here~~
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
function (xmat, y, robdis2 = robdist.hbrfit(xmat), percent = 0.95, 
    intest = ltsreg(xmat, y)$coef) 
{
    xmat = as.matrix(xmat)
    y = as.matrix(y)
    n = dim(xmat)[1]
    p = dim(xmat)[2]
    cut = qchisq(percent, p)
    resids = y - intest[1] - xmat \%*\% as.matrix(intest[2:(p + 
        1)])
    sigma = mad(resids)
    m = psi(cut/robdis2)
    a = resids/(sigma * m)
    c = (median(a) + 3 * mad(a))^2
    h = sqrt(c)/a
    ans = psi(abs(h))
    ans
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
