\name{centerx}
\alias{centerx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Centers a Design Matrix
%%  ~~function to do ... ~~
}
\description{Centers a Design Matrix
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
centerx(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The design matrix which will be centered (all colums sum to 0).
%%     ~~Describe \code{x} here~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
    x = as.matrix(x)
    n = length(x[, 1])
    one = matrix(rep(1, n), ncol = 1)
    x - (one \%*\% t(one)/n) \%*\% x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
