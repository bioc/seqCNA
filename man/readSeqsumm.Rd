\name{readSeqsumm}
\alias{readSeqsumm}
\title{
Function that reads seqsumm summarized files.
}
\description{
The function reads the \code{seqsumm_out.txt} file(s) from the indicated directories and builds a \code{\link{SeqCNAInfo-class}} object, with read count (RC) and genomic information.
}
\usage{
readSeqsumm(build="", tumour.data=NULL, normal.data=NULL, folder=NULL, normal.folder=NULL, resample.win=NULL, sex=TRUE, nproc=2)
}

\arguments{
  \item{build}{
String indicating the genome and build used to generate and annotate the output \code{\link{SeqCNAInfo-class}} object.
Currently, the annotation package supports \code{hg18} and \code{hg19}. This means that common CNV and mappability filters are only available for these builds, and that GC content is estimated from the tumoural sample - or the normal sample if available.
}   
\item{tumour.data}{
A dataframe with the seqsumm information for the tumoural sample.
}
\item{normal.data}{
If applicable, a dataframe with the seqsumm information for the normal sample. Otherwise, disregard this parameter.
}
  \item{folder}{
Path to the folder where the tumoural sample \code{seqsumm_out.txt} file is located.
Only used if no data is passed through the \code{{tumour.data}} parameter.
}
  \item{normal.folder}{
If applicable, path to the folder where the paired normal sample \code{seqsumm_out.txt} file is located. Otherwise, disregard this parameter.
Only used if no data is passed through the \code{{normal.data}} parameter.
}
  \item{resample.win}{
An integer that allows to specify a new bigger summarization window. If used, it must be an exact multiple of the window in the data read.
}
  \item{sex}{
A boolean indicating whether to read sex chromosomes into the output \code{\link{SeqCNAInfo-class}} object.
}
  \item{nproc}{
A value indicating how many processing cores to use for the process of resampling, if applicable. 
Greater values speed up resampling using more CPU cores and RAM memory, but you should not use values greater than the number of cores in your machine. If unsure, the safest value is 1, but most computers nowadays are multi-core, so you could probably go up to 2, 4 or 8.
}
}
\details{
See \code{\link{seqsumm_HCC1143}} for an example table read by the function.
}
\value{
A \code{\link{SeqCNAInfo-class}} object, with information on read count (RC), genome build, and summarization window size and position. If applicable, it also contains paired normal RC. If paired-end mapping (PEM) was used in the alignment, RCs are broken down by read type. 
}
\author{
David Mosen-Ansorena
}
\examples{
data(seqsumm_HCC1143)
rco = readSeqsumm(tumour.data=seqsumm_HCC1143)
}
\keyword{ Summarization }
