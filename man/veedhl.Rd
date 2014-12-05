\name{veedhl}
\alias{veedhl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{veedhl
%%  ~~function to do ... ~~
}
\description{Computes the dispersion-Hodges-Lehmann type variance components.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
veedhl(ehat, center)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ehat}{
%%     ~~Describe \code{ehat} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
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
function (ehat, center) 
{
    numc <- max(center)
    hls <- c(0)
    ehatp <- ehat
    n <- length(ehat)
    for (i in 1:numc) {
        ehatc <- ehat[center == i]
        hls[i] <- hlest(ehatc)
    }
    reff <- hls
    sigb2 <- (pi/3) * ((disp(0, hls, hls, wscores)/numc)^2)
    xmat <- model.matrix(~as.factor(center) - 1)
    ehatp <- ehat - xmat \%*\% hls
    sige2 <- (pi/3) * ((disp(0, ehatp, ehatp, wscores)/n)^2)
    sigt2 <- sigb2 + sige2
    rho <- sigb2/sigt2
    vc <- c(sigt2, sigb2, sige2, rho)
    list(vc = vc, reff = reff)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
