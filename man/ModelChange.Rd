\name{ModelChange}
\alias{joinVCCTest}
\alias{joinECCTest}
\alias{joinCCTest}
\alias{splitVCCTest}
\alias{splitECCTest}
\alias{splitCCTest}
\alias{addECCTest}
\alias{dropECCTest}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model comparison in RCOX models}
\description{
  Methods for altering RCOX models through tests. That is add the most
  significant or drop the most insignificant edge colour class, join the
  two most similar colour classes or split the most heterogeneous
  (composite) colour into atomic colour classes.
}
\usage{
joinVCCTest(m, alpha = 0.05, statistic = c("wald",
                 "dev"), criterion = c("test", "aic", "bic"), details =
                 TRUE)
joinECCTest(m, alpha = 0.05, statistic = c("wald",
                 "dev"), criterion = c("test", "aic", "bic"), details =
                 TRUE)
joinCCTest(m, type = "ecc", alpha = 0.05, statistic = c("wald",
                 "dev"), criterion = c("test", "aic", "bic"), details =
                 TRUE)
splitVCCTest(m,  alpha = 0.05, criterion = c("test",
                 "aic", "bic"), details = TRUE)
splitECCTest(m,  alpha = 0.05, criterion = c("test",
                 "aic", "bic"), details = TRUE)
splitCCTest(m, type = "ecc", alpha = 0.05, criterion = c("test",
                 "aic", "bic"), details = TRUE)
addECCTest(m, alpha = 0.05, criterion = c("test", "aic", "bic"),
                 details = TRUE)
dropECCTest(m, alpha = 0.05, statistic = c("wald", "dev"),
                 criterion = c("test", "aic", "bic"), details = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{An RCOX model object}
  \item{alpha}{Significance level}
  \item{type}{Either 'ecc' (edge colour class) or 'vcc' (vertex colour
    class)}
  \item{details}{If TRUE the details from the parameter comparisons are
    printed. }
  \item{statistic}{Use Wald or likelihood ratio test}
  \item{criterion}{Either 'test', 'aic' or 'bic'}

}
\details{
  splitVCCTest and splitECCTest are simply wrappers for splitCCTest and likewise for
  joinVCCTest, joinECCTest and joinCCTest.
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
bufly0  <- rcox(vcc=list(~me, ~ve, ~al, ~an, ~st), 
             ecc=list(~me:ve, ~me:al, ~ve:al, ~al:an, ~al:st, ~an:st), 
             data=math)

bufly1  <- rcox(~al:an:st, vcc=list(~me+st, ~ve+an), 
                ecc=list(~me:ve+me:al, ~ve:al+al:st),
                data=math)

joinECCTest(bufly0)
splitECCTest(bufly1, criterion='aic')
dropECCTest(bufly1,alpha=0.001,statistic='dev')
addECCTest(bufly1)

}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
