\name{sigma12}
\alias{sigma12}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sigma12
%%  ~~function to do ... ~~
}
\description{Computes square root of matrix Sigma (using Spectral Decomposition).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sigma12(sigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma}{Positive definite matrix.
%%     ~~Describe \code{sigma} here~~
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
function (sigma) 
{
    temp <- eigen(sigma, symmetric = TRUE)
    sigma12 <- temp$vectors \%*\% diag(temp$values^0.5) \%*\% t(temp$vectors)
    sigma12
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
