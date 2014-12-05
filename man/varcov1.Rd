\name{varcov1}
\alias{varcov1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{varcov1
%%  ~~function to do ... ~~
}
\description{Computes asymptotic form  for standard errors as descibed in Section 8.6.4 of Kloke and McKean (2014).

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varcov1(eitbg, center, scores, v12inv, Dmat, wts)
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
  \item{wts}{
%%     ~~Describe \code{wts} here~~
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
function (eitbg, center, scores, v12inv, Dmat, wts) 
{
    numc <- max(center)
    n <- length(eitbg)
    u <- (1:n)/(n + 1)
    sc <- getScores(scores, u)
    rs <- sc[rank(eitbg, ties.method = "first")]
    brrp <- matrix(rep(0, n * n), ncol = n)
    veeinv <- v12inv \%*\% diag(wts) \%*\% v12inv
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
    varcov1 <- solve(part1) \%*\% mid \%*\% solve(part1)
    return(varcov1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
