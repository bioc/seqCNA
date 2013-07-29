\name{SeqCNAInfo-class}
\Rdversion{1.1}
\docType{class}
\alias{SeqCNAInfo-class}
\alias{summary,SeqCNAInfo-method}

\title{Class \code{"SeqCNAInfo"}}
\description{
A structure that contains information on depth of coverage (DOC), genome build, and summarization window size and position. If applicable, it also contains paired normal DOC. If paired-end mapping (PEM) was used in the alignment, DOCs are broken down by read type. Further processing with seqCNA functions of these objects, attach filtering information, and normalized, segmented and called profiles.
}
\section{Objects from the Class}{
Objects can be created by calls to \code{\link{readSeqsumm}} and further modified by \code{\link{applyFilters}}, \code{\link{runSeqnorm}}, \code{\link{runGLAD}} and \code{\link{applyThresholds}}.
}
\section{Slots}{
  \describe{
    \item{\code{tumour}:}{A list with either 1 or 5 numeric vectors, depending on whether the reads are single- or paired-end. Each vector represents the read count for a specific group of read types in a window of size \code{win}. The first vector corresponds to the proper reads.}
    \item{\code{normal}:}{A list analogue to \code{tumour}, but for the normal sample. If there is no information on a normal sample, the list has length zero.}
    \item{\code{seq}:}{A vector indicating the chromosome of each window.}
    \item{\code{pos}:}{A vector indicating the index of each window within the corresponding chromosome.}
    \item{\code{build}:}{Genome and build of the sample. Available builds in the package (hg18, hg19) enable some filters and better GC correction.}
    \item{\code{win}:}{An integer indicating the window size in kilobases.}
    \item{\code{x}:}{The profile against which to perform normalization. Either GC content (estimated or from annotation) or paired-normal read count profile.}
    \item{\code{y}:}{The profile to be normalized, this is, the read count profile of the tumoural sample.}
    \item{\code{skip}:}{A vector with the indexes of the windows that are to be discarded prior to normalization.}
    \item{\code{output}:}{A dataframe with a minimum of 3 columns: chromosome, window start position and normalized profile values. With further processing of the object, columns for the segmented and called profiles are appended.}
    \item{\code{thr}:}{A vector with values that establish the ranges for each copy number.}
    \item{\code{minCN}:}{The copy number of the lowest range defined by the thresholds.}
    \item{\code{gc}:}{Mean GC content of the reads within each genomic window. The means from the tumoural sample are used unless there is a paired-normal sample.}
    \item{\code{mapq}:}{Mean mapping quality of the reads within each genomic window. The means from the tumoural sample are used unless there is a paired-normal sample.}
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "SeqCNAInfo")}: Prints on screen a summary of the object and the applied processing.}
	 }
}
\author{
David Mosen-Ansorena
}
\examples{
showClass("SeqCNAInfo")
}
\keyword{classes}
