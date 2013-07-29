\name{applyFilters}
\alias{applyFilters}
\title{
The function applies a set of filters to the raw profile.
}
\description{
The available filters are: PEM-based (only for paired-end data), trimming of extreme values, mappability-based and common CNVs.
}
\usage{
applyFilters(rco, pem.filter=0, trim.filter=0, mapp.filter=0, mapq.filter=0, cnv.filter=FALSE, plots=TRUE, folder=NULL, nproc=2)
}
\arguments{
  \item{rco}{
A \code{\link{SeqCNAInfo-class}} object, with read count (RC) and genomic information, normally the output of the \code{\link{readSeqsumm}} function.
}
  \item{pem.filter}{
seqsumm summarizes reads not only by window, but also by read type (i.e. SAM read flags). The ratio between proper and improper reads is calculated for each window, where higher ratios correlate with windows in centromeres, structural polymorphisms,... and in general regions advisable to be filtered.
The value indicates the quantile of windows to filter based on the improper read ratio.
}
  \item{trim.filter}{
A numeric value or a vector with two numeric values. Set to a value between 0 and 1 to filter such top and bottom quantiles (e.g. 0.01). If two numeric values are provided, the first one is used for the upper quantile and the second for the lower quantile.
If the normal sample is not provided, set to 1 for automatic upper threshold selection (lower thresohold can be set with a second value or take the same value as the upper threshold).
Pre-normalization is performed over the paired-normal read count, or the tumoural one if the paired-normal is not available, against the GC content. This keeps  the trimming from being biased by the influence of GC content on the read count.
}
  \item{mapp.filter}{
Only available if the annotation package has information on the sample genome and build (currently hg18 and hg19). Mappability ranges from 0 to 1. The higher the parameter value, the less windows will be filtered.
}
  \item{mapq.filter}{
This filter discards genomic windows with low mean mapping quality in the reads of the tumoural sample. The mean mapping quality in a window may range from 0 to a few tens, with an upper limit of 255. Therefore, the higher the parameter value, the more windows will be filtered.
}
  \item{cnv.filter}{
A boolean indicating whether to apply the CNV-based filter. If applied, windows with a common CNV (Altshuler et al., 2010) spanning at least 95 per cent of the window are filtered.
This filter is available for completion purposes, but, given that the frequency of common CNVs can be as low as 0.01, very few polymorphisms are captured in comparison to the amount of expected false positives.
}
\item{plots}{
A boolean indicating whether to output plots to the \code{folder} folder.
}
  \item{folder}{
If the \code{plots} parameter is set to TRUE, path to the folder where the plots with filtering information are to be generated.
If no folder is indicated or does not exist, plots will be displayed within R.
}
  \item{nproc}{
A value indicating how many processing cores to use for the processes of building genome information and automatic trimming, if applicable. 
Greater values speed up the processsing using more CPU cores and RAM memory, but you should not use values greater than the number of cores in your machine. If unsure, the safest value is 1, but most computers nowadays are multi-core, so you could probably go up to 2, 4 or 8.
}
}
\details{
The function outputs two plots to the \code{folder} folder. We suggest reviewing these plots in order to adjust the filters to adequate values.
For the Venn diagrams, functions by Thomas Girke are used (see http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html).
}
\value{
A \code{\link{SeqCNAInfo-class}} object, with additional information on which windows to filter and the profiles needed for subsequent normalization in \code{\link{runSeqnorm}}.
}
\references{
Integrating common and rare genetic variation in diverse human populations.
Altshuler DM, Gibbs RA, Brooks LD, McEwen JE.
Nature. 2010 Sep 2; 467:52-8
}
\author{
David Mosen-Ansorena
}
\examples{
data(seqsumm_HCC1143)
rco = readSeqsumm(tumour.data=seqsumm_HCC1143)

### FILTERING ###

pem.filter = 0
trim.filter = 1
cnv.filter = FALSE
mapp.filter = 0
mapq.filter = 2

rco = applyFilters(rco, pem.filter, trim.filter, mapp.filter, mapq.filter, cnv.filter, plots=FALSE)
}

\keyword{ Filtering }
