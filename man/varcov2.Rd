\name{varcov2}
\alias{varcov2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{varcov2
%%  ~~function to do ... ~~
}
\description{Computes the approximate type estimator for standard errors as descibed in Section 8.6.4 of Kloke and McKean (2014).
This is the default method.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varcov2(eitbg, center, scores, v12inv, Dmat, p, delta, hparm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eitbg}{
%%     ~~Describe \code{eitbg} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
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
function (eitbg, center, scores, v12inv, Dmat, p, delta, hparm) 
{
    numc <- max(center)
    n <- length(eitbg)
    u <- (1:n)/(n + 1)
    sc <- getScores(scores, u)
    rs <- sc[rank(eitbg, ties.method = "first")]
    brrp <- matrix(rep(0, n * n), ncol = n)
    tauhat <- gettauF0(eitbg, p, scores, delta, hparm)
    veeinv <- v12inv \%*\% v12inv
    part1 <- t(Dmat) \%*\% veeinv \%*\% Dmat
    i1 <- 1
    i2 <- 1
    for (i in 1:numc) {
        rsc <- rs[center == i]
        ni <- length(rsc)
        rscbar <- mean(rsc)
        rsc <- rsc - rscbar
        rrp <- rsc \%*\% t(rsc)
        i2 <- i1 + ni - 1
        brrp[i1:i2, i1:i2] <- rrp
        i1 <- i1 + ni
    }
    mid <- t(Dmat) \%*\% v12inv \%*\% brrp \%*\% v12inv \%*\% Dmat
    varcov2 <- (tauhat^2) * solve(part1) \%*\% mid \%*\% solve(part1)
    return(varcov2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
