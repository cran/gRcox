\name{tr}
\alias{trX}
\alias{trXY}
\alias{trXYZ}
\alias{trXYXY}
\alias{trXYZY}

\title{Calculate trace of various matrix products}
\description{
  Calculate trace of various matrix products.
}
\usage{
trX(x, PACKAGE = "gRcox")
trXY(x, y, PACKAGE = "gRcox")
trXYXY(x, y, PACKAGE = "gRcox")
trXYZ(x, y, z, PACKAGE = "gRcox")
trXYZY(x, y, z, PACKAGE = "gRcox")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y,z}{Compatible matrices }
  \item{PACKAGE}{"gRcox"}
}
%\details{
%
%}
\value{
  A number
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
\note{
  The storage mode of the matrices must be 'double', see examples below
}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

}
\keyword{utilities}% at least one, from doc/KEYWORDS
