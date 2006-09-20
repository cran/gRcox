\name{display}
\alias{display}
\alias{display.rcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays an RCOX model graphically (using the Rgraphviz package)}
\description{
  The display function creates a coloured graph representing an RCOX model.
}
\usage{
display(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An RCOX model}
}
\details{
  Note that the Rgraphviz package must be loaded.
}
\value{
  NULL
}
%\references{.....}
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{ ~~further notes~~ 
%
%% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
\dontrun{
data(math)
bufly1  <- rcox(~al:an:st, vcc=list(~me+st, ~ve+an), 
                ecc=list(~me:ve+me:al, ~ve:al+al:st),
                data=math)
library(Rgraphviz)
display(bufly1)
}
}
\keyword{hplot}% at least one, from doc/KEYWORDS

