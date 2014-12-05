\name{getAp}
\alias{getAp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getAp
%%  ~~function to do ... ~~
}
\description{If the link is the identity, i.e., a linear model then the function aprime is
the design matrix times beta.  If another link is to be used hen this function needs to be changed.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getAp(xmat, center, beta, geemod = "LM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmat}{Design matrix.
%%     ~~Describe \code{xmat} here~~
}
  \item{center}{Center (block) indicator.
%%     ~~Describe \code{center} here~~
}
  \item{beta}{Current value of beta
%%     ~~Describe \code{beta} here~~
}
  \item{geemod}{Currently only geemod="LM" functions.
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
        apbeta <- xmat \%*\% beta
    }
    return(apbeta)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
