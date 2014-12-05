\name{hlest}
\alias{hlest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hodges-Lehmann (median of Walsh averages) Estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hlest(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Sample
%%     ~~Describe \code{x} here~~
}
}
\details{If med is set to FALSE in geerfit call then the Hodges-Lehmann (estimator is computed.
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
\author{Joseph McKean
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
function (x) 
{
    mat <- pairup(x, type = "leq")
    hlest <- median((mat[, 1] + mat[, 2])/2)
    return(hlest)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
