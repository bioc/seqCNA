\name{writeCNProfile}
\alias{writeCNProfile}
%- Also NEED an '\alias' for EACH other topic rcoumented here.
\title{
Outputs the processed profile of a \code{\link{SeqCNAInfo-class}} object to the specified folder.
}
\description{
The output columns depend on the processing applied on the \code{\link{SeqCNAInfo-class}} object.
}
\usage{
writeCNProfile(rco, folder)
}
\arguments{
  \item{rco}{
A \code{\link{SeqCNAInfo-class}} object.
}
  \item{folder}{
Path to the folder where the table with the output from the \code{\link{SeqCNAInfo-class}} object is to be generated.
}
}
\value{
Nothing is returned from this function. Check the \code{folder} folder for a text file called \code{seqnorm_out.txt}.
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
rco = applyThresholds(rco, seq(-0.9,4,by=0.9), 1)

writeCNProfile(rco, tempdir())
}
\keyword{ Output }
