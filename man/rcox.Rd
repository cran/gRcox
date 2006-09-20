\name{rcox}
\alias{rcox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Main function for specifying RCON/RCOR models}
\description{
  This is the main function for specifying and fitting RCON/RCOR models in the
  package. 
}
\usage{
rcox(gm = NULL, vcc = NULL, ecc = NULL, type = c("rcon", "rcor"),
method = c("scoring", "ipm", "hyd", "user"),
fit = TRUE,
data = NULL, S = NULL, n = NULL, Kstart,
control = rcox.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gm}{Generating class for a grapical Gaussian model, see
    'Examples' for an illustration}
  \item{vcc}{List of vertex colour classes for the model}
  \item{ecc}{List of edge colour classes for the model}
  \item{type}{Type of model. Default is RCON}
  \item{method}{Estimation method. Default is 'scoring' which is
    stabilised Fisher scoring. An alternative is 'ipm' which is
    iterative partial maximisation. The methods 'hyd' and 'user' are for
  internal use and should not be called directly}
  \item{fit}{Should the model be fitted}
  \item{data}{A dataframe}
  \item{S}{An empirical covariance matrix (as alternative to giving data
  as a dataframe)}
  \item{n}{The number of observations (which is needed if data is
    specified as an empirical covariance matrix)}
  \item{Kstart}{An initial value for K}
  \item{control}{Controlling the fitting algorithms}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  A model object of type 'rcox'.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{ ~~further notes~~ 

%% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

}
\keyword{models}% at least one, from doc/KEYWORDS
