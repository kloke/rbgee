\name{varcov3}
\alias{varcov3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{varcov3
%%  ~~function to do ... ~~
}
\description{Computes a sandwich type estimator for standard errors as descibed in Section 8.6.4 of Kloke and McKean (2014).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varcov3(eitbg, scores, v12inv, Dmat, p, delta, hparm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eitbg}{
%%     ~~Describe \code{eitbg} here~~
}
  \item{scores}{
%%     ~~Describe \code{scores} here~~
}
  \item{v12inv}{
%%     ~~Describe \code{v12inv} here~~
}
  \item{Dmat}{
%%     ~~Describe \code{Dmat} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{hparm}{
%%     ~~Describe \code{hparm} here~~
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
function (eitbg, scores, v12inv, Dmat, p, delta, hparm) 
{
    n <- length(eitbg)
    u <- (1:n)/(n + 1)
    sc <- getScores(scores, u)
    tauhat <- gettauF0(eitbg, p, scores, delta, hparm)
    veeinv <- v12inv \%*\% v12inv
    part1 <- t(Dmat) \%*\% veeinv \%*\% Dmat
    varcov3 <- (tauhat^2) * solve(part1)
    return(varcov3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
