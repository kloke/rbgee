\name{wtmat}
\alias{wtmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{wtmat
%%  ~~function to do ... ~~
}
\description{Computes the weight matrix W used in the rank-based GEE fit.
See Section 8.6 of Kloke and McKean (2014).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wtmat(dmatb, eitb, med = TRUE, scores = wscores, hbrs = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmatb}{
%%     ~~Describe \code{dmatb} here~~
}
  \item{eitb}{
%%     ~~Describe \code{eitb} here~~
}
  \item{med}{
%%     ~~Describe \code{med} here~~
}
  \item{scores}{
%%     ~~Describe \code{scores} here~~
}
  \item{hbrs}{
%%     ~~Describe \code{hbrs} here~~
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
function (dmatb, eitb, med = TRUE, scores = wscores, hbrs = FALSE) 
{
    n <- length(eitb)
    u <- (1:n)/(n + 1)
    eps <- 1e-06
    sc <- getScores(scores, u)
    top <- sc[rank(eitb, ties.method = "first")]
    wts <- rep(0, n)
    if (med) {
        mb <- median(eitb)
        ehat <- eitb - mb
        wts <- top/ehat
        wts[is.infinite(wts)] <- 0
        wts[is.nan(wts)] <- 0
        wts[wts == 0] <- max(wts)
    }
    else {
        scl <- sc[sc < 0]
        mp1 <- length(scl) + 1
        eitbs <- sort(eitb)
        mb <- eitbs[mp1]
        ehat <- eitb - mb
        wts <- top/ehat
        wts[is.infinite(wts)] <- 0
        wts[is.nan(wts)] <- 0
        wts[wts == 0] <- max(wts)
    }
    if (hbrs) {
        wts2 <- hbrwtsnew(dmatb, eitb)
        newwts <- wts * wts2
        newwts[is.infinite(wts)] <- 0
        newwts[is.nan(wts)] <- 0
        newwts[wts == 0] <- max(wts)
        wts <- newwts
    }
    list(wts = wts, mb = mb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
