\name{refGenchr1}
\alias{refGenchr1}
\docType{data}
\title{
Gene annotation of chromosome 1 of human reference genome.
}
\description{
This data is a part of the gene annotation of chr1 of human reference genome. 
}
\usage{data("refGenchr1")}
\format{
  A data frame with 3066 observations on the following 9 variables.
  \describe{
    \item{\code{V1}}{Name of chromosome}
    \item{\code{V2}}{Source}
    \item{\code{V3}}{Feature}
    \item{\code{V4}}{Start position}
    \item{\code{V5}}{End position}
    \item{\code{V6}}{Score}
    \item{\code{V7}}{Strand}
    \item{\code{V8}}{Frame}
    \item{\code{V9}}{Attribute}
  }
}
\details{
This is an example annotation file of chr1 (human reference genome) whose start position 
is less than or equal 2000000.
}
\source{
Downloaded from UCSC genome browser.
}
\references{
UCSC genome browser.
}
\examples{
data(refGenchr1)
## maybe str(refGenchr1) ; plot(refGenchr1) ...
}
\keyword{datasets}
