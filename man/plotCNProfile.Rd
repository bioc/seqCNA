\name{plotCNProfile}
\alias{plotCNProfile}
\title{
Plots the profile of a \code{\link{SeqCNAInfo-class}} object.
}
\description{
The plotted elements depend on the processing applied on the \code{\link{SeqCNAInfo-class}} object.
}
\usage{
plotCNProfile(rco, folder = NULL)
}
\arguments{
  \item{rco}{
A \code{\link{SeqCNAInfo-class}} object.
}
  \item{folder}{
Path to the folder where the plots with the output from the \code{\link{SeqCNAInfo-class}} object are to be generated.
If no folder is indicated or does not exist, plots will be displayed within R.
}
}
\value{
Nothing is returned from this function. Check the \code{folder} folder for a file called \code{seqnorm_out.jpg}.
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

plotCNProfile(rco)
}
\keyword{ Output }
