\name{veemm}
\alias{veemm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{veemm
%%  ~~function to do ... ~~
}
\description{Computes the mad-median type of variance component estimator.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
veemm(ehat, center)
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
    meds <- c(0)
    ehatp <- ehat
    n <- length(ehat)
    for (i in 1:numc) {
        ehatc <- ehat[center == i]
        meds[i] <- median(ehatc)
    }
    reff <- meds
    sigb2 <- mad(meds)^2
    xmat <- model.matrix(~as.factor(center) - 1)
    ehatp <- ehat - xmat \%*\% meds
    sige2 <- mad(ehatp)^2
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
