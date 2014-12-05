\name{getD}
\alias{getD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getD
%%  ~~function to do ... ~~
}
\description{Computes the Hessian matrix D; see Section 8.6 of Kloke and McKean (2014).
Currently the only option is the identity link; i.e., D is the design matrix.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getD(xmat, center, beta, geemod = "LM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmat}{
%%     ~~Describe \code{xmat} here~~
}
  \item{center}{
%%     ~~Describe \code{center} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{geemod}{
%%     ~~Describe \code{geemod} here~~
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
function (xmat, center, beta, geemod = "LM") 
{
    if (geemod == "LM") {
        D <- xmat
    }
    return(D)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
