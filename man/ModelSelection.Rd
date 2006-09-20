\name{ModelSelection}
\alias{stepJoinVCC}
\alias{stepJoinECC}
\alias{stepJoinCC}
\alias{stepSplitVCC}
\alias{stepSplitECC}
\alias{stepSplitCC}
\alias{stepAddECC}
\alias{stepDropECC}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stepwise selection of RCOX models}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stepJoinVCC(m,  alpha = 0.05, details=TRUE,
                 criterion = c("test", "aic", "bic"), statistic =
                 c("wald", "dev"), display = FALSE)
stepJoinECC(m,  alpha = 0.05, details=TRUE,
                 criterion = c("test", "aic", "bic"), statistic =
                 c("wald", "dev"), display = FALSE)
stepJoinCC(m, type = "ecc", alpha = 0.05, details=TRUE,
                 criterion = c("test", "aic", "bic"), statistic =
                 c("wald", "dev"), display = FALSE)
stepSplitVCC(m,  alpha = 0.05,details=TRUE,
                 criterion = c("test", "aic", "bic"), display = FALSE)
stepSplitECC(m,  alpha = 0.05,details=TRUE,
                 criterion = c("test", "aic", "bic"), display = FALSE)
stepSplitCC(m, type = "ecc", alpha = 0.05,details=TRUE,
                 criterion = c("test", "aic", "bic"), display = FALSE)
stepDropECC(m, alpha = 0.05, statistic = c("wald", "dev"),
                 criterion = c("test", "aic", "bic"), details = TRUE,  display = FALSE)
stepAddECC(m, alpha = 0.05, criterion = c("test", "aic", "bic"), details = TRUE, display = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{An RCOX model object}
  \item{alpha}{Significance level}
  \item{type}{Either 'ecc' (edge colour class) or 'vcc' (vertex colour class)}
  \item{details}{Show details of the tests}
  \item{statistic}{Use Wald or likelihood ratio test}
   \item{criterion}{Either 'test', 'aic' or 'bic'}
   \item{display}{If TRUE and if the Rgraphviz package is loaded then
     the models are displayed graphically.}
 }
\details{
  stepSplitVCC and stepSplitECC are simply wrappers for stepSplitCC and likewise for
  stepJoinVCC, stepJoinECC and stepJoinCC.

}
\value{
  An RCOX model object.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Søren Højsgaard, sorenh@agrsci.dk}
%\note{ ~~further notes~~ %
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(math)
bufly0  <- rcox(vcc=list(~me, ~ve, ~al, ~an, ~st), 
             ecc=list(~me:ve, ~me:al, ~ve:al, ~al:an, ~al:st, ~an:st), 
             data=math)

bufly01<-stepJoinVCC(bufly0,statistic='dev')

bufly02<-stepJoinECC(bufly01, details=FALSE, statistic='dev')

bufly2  <- rcox(vcc=list(~al, ~me+st, ~ve+an), 
              ecc=list(~me:ve+me:al+ve:al, ~al:an+al:st+an:st), 
              data=math)
stepSplitECC(bufly2)
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
