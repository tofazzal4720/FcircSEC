\name{chr1}
\alias{chr1}
\docType{data}
\title{
chr1 of Reference genome 
}
\description{
part of chr1 (first 2000000 nucleotides) from the reference genome. 
}
\usage{data("chr1")}
\format{
  A data frame with 1 observations on the following 2 variables.
  \describe{
    \item{\code{seq_name}}{Name of the chromosome}
    \item{\code{sequence}}{Nucleotides sequence}
  }
}
\details{
This data contains only the first 2000000 nucleotides of the reference 
chromosome 1.  
}
\source{
The full reference genome was downloaded from UCSC.
}
\references{
UCSC genome browser.
}
\examples{
data(chr1)
## maybe str(chr1) ; plot(chr1) ...
}
\keyword{datasets}
