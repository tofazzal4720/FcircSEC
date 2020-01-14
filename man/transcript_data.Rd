\name{transcript_data}
\alias{transcript_data}
\docType{data}
\title{
Transcript data obtained from gene annotation file.
}
\description{
This data is produced by the r package FcircSEC.
}
\usage{data("transcript_data")}
\format{
  A data frame with 177 observations on the following 9 variables.
  \describe{
    \item{\code{transcript_id}}{Transcript name}
	\item{\code{chr}}{Chromosome}
    \item{\code{strand}}{Strand}
    \item{\code{trans_start}}{Transcript start}
    \item{\code{trans_end}}{Transcript end}
    \item{\code{exon_count}}{Number of exons}
    \item{\code{exon_starts}}{Start position of exons}
	\item{\code{exon_ends}}{End position of exons}
	\item{\code{gene}}{Gene name}
  }
}
\details{
This is an example data produced by the function transcriptExtract of the r package FcircSEC.
}
\source{
This data is produced by r package FcircSEC.
}
\references{
R package FcircSEC.
}
\examples{
data(transcript_data)
## maybe str(transcript_data) ; plot(transcript_data) ...
}
\keyword{datasets}
