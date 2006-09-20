\name{Model Editing}
\alias{joinVCC}
\alias{joinECC}
\alias{joinCC}
\alias{splitVCC}
\alias{splitECC}
\alias{splitCC}
\alias{addECC}
\alias{dropECC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Editing of RCOX models}
\description{
  Functions for editing (i.e. modifying) RCOX models.
}
\usage{
joinVCC(cc, m,  fit = TRUE)
joinECC(cc, m,  fit = TRUE)
joinCC(cc, m, type = "ecc", fit = TRUE)
splitVCC(cc, m, fit = TRUE)
splitECC(cc, m, fit = TRUE)
splitCC(cc, m, type = "ecc", fit = TRUE) 
addECC(cc, m, fit = TRUE)
dropECC(cc, m, fit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cc}{A colour class}
  \item{m}{An RCOX model}
  \item{type}{Either 'vcc' for vertex colour class or 'ecc' for edge
    colour class}
  \item{fit}{If FALSE the new model is not fitted.}
}
\details{
  splitVCC and splitECC are simply wrappers for splitCC and likewise for
  joinVCC, joinECC and joinCC.
}
\value{
  An RCOX model object.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(math)
bufly1  <- rcox(~al:an:st, vcc=list(~me+st, ~ve+an), 
                ecc=list(~me:ve+me:al, ~ve:al+al:st),
                data=math)

joinVCC(list(~me+st, ~ve+an), bufly1)
joinECC(list(~al:an, ~an:st), bufly1)
splitCC(~ve+an,bufly1,type='vcc')
splitCC(~me:ve+me:al,bufly1,type='ecc')
splitVCC(~ve+an, bufly1)
splitECC(~me:ve+me:al, bufly1)
addECC(~me:an+ve:st,bufly1)
dropECC(~me:ve+me:al,bufly1)
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
