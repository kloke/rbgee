\name{veemat}
\alias{veemat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{veemat
%%  ~~function to do ... ~~
}
\description{Computes the working covariance matrix.   The method is set by the arguments structure and substructure.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
veemat(ehat, center, structure = "CS", substructure = "MM", scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ehat}{
%%     ~~Describe \code{ehat} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
}
  \item{structure}{"WI" for working independence; "CS" for compound symmetry (default); "AR" for autoregressive 1.
%%     ~~Describe \code{structure} here~~
}
  \item{substructure}{"MM" for mad-median option of variance components (default); "DHL" for disoersion-HL option.
%%     ~~Describe \code{substructure} here~~
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
function (ehat, center, structure = "CS", substructure = "MM", 
    scores) 
{
    numc <- max(center)
    vc <- c()
    n <- length(ehat)
    v12 <- matrix(rep(0, n^2), ncol = n)
    v12inv <- matrix(rep(0, n^2), ncol = n)
    if (structure == "CS") {
        if (substructure == "DHL") {
            vc <- veedhl(ehat, center)$vc
        }
        if (substructure == "MM") {
            vc <- veemm(ehat, center)$vc
        }
        i1 <- 1
        i2 <- 1
        for (i in 1:numc) {
            esee = ehat[center == i]
            ni <- length(esee)
            i2 <- i1 + ni - 1
            vt <- vc[1] * ((1 - vc[4]) * diag(rep(1, ni)) + vc[4] * 
                rep(1, ni) \%*\% t(rep(1, ni)))
            vt12 <- sigma12(vt)
            vt12inv <- sigma12inv(vt)
            v12[i1:i2, i1:i2] <- vt12
            v12inv[i1:i2, i1:i2] <- vt12inv
            i1 <- i1 + ni
        }
    }
    if (structure == "WI") {
        v12 <- diag(rep(1, n))
        v12inv <- diag(rep(1, n))
    }
    if (structure == "AR") {
        vc <- veear1m(ehat, center, scores)
        i1 <- 1
        i2 <- 1
        for (i in 1:numc) {
            esee = ehat[center == i]
            ni <- length(esee)
            i2 <- i1 + ni - 1
            kap <- vc[1]/(1 - vc[2]^2)
            vt <- matrix(rep(0, ni^2), ncol = ni)
            for (ip in 1:ni) {
                for (jp in 1:ni) {
                  vt[ip, jp] <- kap * vc[2]^abs((jp - ip))
                }
            }
            vt12 <- sigma12(vt)
            vt12inv <- sigma12inv(vt)
            v12[i1:i2, i1:i2] <- vt12
            v12inv[i1:i2, i1:i2] <- vt12inv
            i1 <- i1 + ni
        }
    }
    list(v12 = v12, v12inv = v12inv, vc = vc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
