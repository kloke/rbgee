\name{gencs2}
\alias{gencs2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gencs2
%%  ~~function to do ... ~~
}
\description{Simulator for a GEE model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gencs2(vecni, vc, beta, xmat, itypee = 1, itypeb = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vecni}{
%%     ~~Describe \code{vecni} here~~
}
  \item{vc}{
%%     ~~Describe \code{vc} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{xmat}{
%%     ~~Describe \code{xmat} here~~
}
  \item{itypee}{
%%     ~~Describe \code{itypee} here~~
}
  \item{itypeb}{
%%     ~~Describe \code{itypeb} here~~
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
function (vecni, vc, beta, xmat, itypee = 1, itypeb = 1) 
{
    numc <- length(vecni)
    n <- sum(vecni)
    center <- rep(0, n)
    ind <- rep(0, n)
    errs <- rep(0, n)
    i1 <- 1
    i2 <- 0
    for (j in 1:numc) {
        if (itypeb == 1) {
            bj <- rnorm(1, 0, sqrt(vc[1]))
        }
        i2 <- i2 + vecni[j]
        ni <- vecni[j]
        ni2 <- floor(ni/2)
        center[i1:i2] <- rep(j, ni)
        if (itypee == 1) {
            ei <- rnorm(ni, 0, sqrt(vc[2]))
        }
        errs[i1:i2] <- ei + bj
        i1 <- i1 + ni
    }
    y <- xmat \%*\% beta + errs
    list(y = y, center = center)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
