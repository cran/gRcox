\name{Parameter Comparison}
\alias{compareVCC}
\alias{compareECC}
\alias{compareCC}
\alias{compareAllVCC}
\alias{compareAllECC}
\alias{compareAllCC}
\alias{zeroECC}
\alias{zeroAllECC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise comparison of parameters in RCOX models}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compareVCC(cc1, cc2, mf, statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
compareECC(cc1, cc2, mf, statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
compareCC(cc1, cc2, mf, type = "ecc", statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
compareAllVCC(mf, order = TRUE, statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
compareAllECC(mf, order = TRUE, statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
compareAllCC(mf, type = "ecc", order = TRUE, statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
zeroECC(ecc, mf, statistic = c("wald", "dev"), criterion =  c("test", "aic", "bic"))
zeroAllECC (mf, order = TRUE, statistic = c("wald", "dev"), criterion = c("test", "aic", "bic"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cc1}{A colour class}
  \item{cc2}{A colour class}
  \item{ecc}{An edge colour class}
  \item{mf}{An RCOX model object}
  \item{type}{Either 'ecc' (edge colour class) or 'vcc' (vertex colour class)}
  \item{order}{Should the result be sorted.}
    \item{statistic}{Use Wald or likelihood ratio test}
  \item{criterion}{Either 'test', 'aic' or 'bic'}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A dataframe with the test result.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
