\name{applyThresholds}
\alias{applyThresholds}
\title{
Apply thresholds on the segmented profile in order to call copy numbers.
}
\description{
Given a set of thresholds and the copy number of the lower resulting range, copy numbers are called.
}
\usage{
applyThresholds(rco, thresholds, min.CN)
}
\arguments{
  \item{rco}{
A \code{\link{SeqCNAInfo-class}} object, with read count (RC) and genomic information, normally the output of the \code{\link{runGLAD}} function.
}
  \item{thresholds}{
A vector with values that establish the ranges for each copy number.
}
  \item{min.CN}{
The copy number of the lowest range defined by the thresholds.
}
}
\value{
A \code{\link{SeqCNAInfo-class}} object, with additional information on the copy numbers called for each genomic window.
}
\author{
David Mosen-Ansorena
}
\examples{
data(seqsumm_HCC1143)
rco = readSeqsumm(tumour.data=seqsumm_HCC1143)
rco = applyFilters(rco, 0, 1, 0, 2, FALSE, plots=FALSE)
rco = runSeqnorm(rco, plots=FALSE)
rco = runGLAD(rco)

### CALLING ###

thresholds = seq(-0.9,4,by=0.9)
min.CN = 1

rco = applyThresholds(rco, thresholds, min.CN)
}
\keyword{ Calling }
