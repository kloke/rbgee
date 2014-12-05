\name{bigstars}
\alias{bigstars}
\docType{data}
\title{Jittered Stars Data
%%   ~~ data name/kind ... ~~
}
\description{This is an expanded version of the famous stars data; see, for example,
Hettmansperger and McKean (2011).
The first 47 data points are the original stars data.
A second copy of the x values,  a jittered version of the first
47 y values, and a block variable (1 if in first 47, else 2) constitute the data set.
It is used to check high breakdown GEE fits.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(bigstars)}
\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{bigx}}{a numeric vector}
    \item{\code{bigy}}{a numeric vector}
    \item{\code{block}}{a numeric vector}
  }
}
\details{The first 47 x's are log temperature of stars
while the first 47 y values are the log intensity of light.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(bigstars)
## maybe str(bigstars) ; plot(bigstars) ...
}
\keyword{datasets}
