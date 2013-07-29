\name{runGLAD}
\alias{runGLAD}
\title{
Wrapper function that runs the GLAD segmentation algorithm on a \code{\link{SeqCNAInfo-class}} object.
}
\description{
The function calls the corresponding function in the \code{GLAD} library with appropriate parameters.
}
\usage{
runGLAD(rco, lambdabreak=8, nproc=2)
}
\arguments{
  \item{rco}{
A \code{\link{SeqCNAInfo-class}} object, with read count (RC) and genomic information, normally the output of the \code{\link{runSeqnorm}} function.
}
  \item{lambdabreak}{
Penalty term used in GLAD during the optimization of the number of breakpoints step. Higher values result in lower amount of segments.
}
  \item{nproc}{
A value indicating how many processing cores to use in the segmentation. 
Greater values speed up segmentation using more CPU cores and RAM memory, but you should not use values greater than the number of cores in your machine. If unsure, the safest value is 1, but most computers nowadays are multi-core, so you could probably go up to 2, 4 or 8.
}
}
\value{
A \code{\link{SeqCNAInfo-class}} object, with additional information on the segmented (i.e. smoothed) profile.
}
\references{
Analysis of array CGH data: from signal ratio to gain and loss of DNA regions.
Hupe P, Stransky N, Thiery JP, Radvanyi F, Barillot E.
Bioinformatics. 2004 Dec 12; 20(18):3413-22.
}
\author{
David Mosen-Ansorena
}
\examples{
data(seqsumm_HCC1143)
rco = readSeqsumm(tumour.data=seqsumm_HCC1143)
rco = applyFilters(rco, 0, 1, 0, 2, FALSE, plots=FALSE)
rco = runSeqnorm(rco, rco@win, "quadratic", nproc=2, plots=FALSE)

### SEGMENTATION ###

rco = runGLAD(rco)
}
\keyword{ Segmentation }
