\name{runSeqnorm}
\alias{runSeqnorm}
\title{
This function calls seqnorm, which normalizes the tumoural profile with respect to either a paired sample or GC content.
}
\description{
The seqnorm method performs several regressions on homogeneous regions of the tumoural genomic profile.
}
\usage{
runSeqnorm(rco, norm.win = NULL, method = "quadratic", nproc = 2, plots = TRUE, folder = NULL)
}
\arguments{
  \item{rco}{
A \code{\link{SeqCNAInfo-class}} object, with read count (RC) and genomic information, normally the output of the \code{\link{readSeqsumm}} function.
}
  \item{norm.win}{
The genomic window size, given in kilobases, in which to normalize the read count profile.
Leave it to the default NULL value for automatic selection.
The window size should be the same as or a multiple of the summarization window in \code{\link{runSeqsumm}}. If not a multiple, the window size will be adjusted to the nearest one.
For efficiency, we recommend to use the same value as the summarization window size, unless it is very small (e.g. 2Kbp).
}
  \item{method}{
A string being \code{"loess"}, \code{"cubic"} or \code{"quadratic"}, indicating which form of regression to use in the normalization.
LOESS is more sensitive to noise and slower than the other two, which apply polynomial regression. In principle, we recommend using a quadratic polynomial regression, but there might be ocassions in which the others yield better fits.
}
  \item{nproc}{
A value indicating how many processing cores to use in the segmentation prior to regression and in the segment regressions. 
Greater values speed up normalization using more CPU cores and RAM memory, but you should not use values greater than the number of cores in your machine. If unsure, the safest value is 1, but most computers nowadays are multi-core, so you could probably go up to 2, 4 or 8.
}
  \item{plots}{
A boolean value indicating whether to output normalization plots to the \code{folder}. These are useful for checking the correctness of the normalization and the improvements over the standard normalization.
}
  \item{folder}{
If the \code{plots} parameter is set to TRUE, path to the folder where the plots with normalization information are to be generated.
If no folder is indicated or does not exist, plots will be displayed within R.
}
}
\value{
A \code{\link{SeqCNAInfo-class}} object, with additional information on the normalized profile.
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

### NORMALIZATION ###

rco = runSeqnorm(rco, plots=FALSE)
}
\keyword{ Normalization }
