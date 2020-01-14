\name{output_CIRI}
\alias{output_CIRI}
\docType{data}
\title{
Output obtained from circRNA prediction tool.
}
\description{
Four columns of the circRNA prediction tool.
}
\usage{data("output_CIRI")}
\format{
  A data frame with 14 observations on the following 4 variables.
  \describe{
    \item{\code{V1}}{Chromosome name}
    \item{\code{V2}}{Start position}
    \item{\code{V3}}{End position}
    \item{\code{V4}}{Strand}
  }
}
\details{
Four columns of the circRNA prediction tool. The columns
have information about chromosome, circRNA start and end positions, and strand.
}
\source{
An example data from the circRNA prediction tool CIRI.
}
\references{
This data is obtained from the circRNA prediction tool CIRI.
}
\examples{
data(output_CIRI)
## maybe str(output_CIRI) ; plot(output_CIRI) ...
}
\keyword{datasets}
