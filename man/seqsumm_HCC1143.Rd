\name{seqsumm_HCC1143}
\alias{seqsumm_HCC1143}
\docType{data}
\title{
Example of a seqsumm output, namely for a sequenced sample of the HCC1143 cancer cell-line.
}
\description{
The \code{\link{runSeqsumm}} function was run over a SAM file at a summarization window of 200Kbps and only chromosomes 1 to 5 were kept.
}
\usage{data(seqsumm_HCC1143)}
\format{
  A data frame with 5314 genomic windows on the following 5 variables.
  \describe{
    \item{\code{chrom}}{A numeric vector with the chromosome of the corresponding window.}
    \item{\code{win.start}}{A numeric vector with the first zero-based base pair index of the corresponding window.}
    \item{\code{reads.gc}}{A numeric vector with the mean GC content for the reads within the corresponding window.}
    \item{\code{reads.mapq}}{A numeric vector with the mean mapping quality for the reads within the corresponding window.}
    \item{\code{counts}}{A numeric vector with the read count within the corresponding window.}
  }
}
\references{
High-resolution mapping of copy-number alterations with massively parallel sequencing.
Chiang DY, Getz G, Jaffe DB, O'Kelly MJ, Zhao X, Carter SL, Russ C, Nusbaum C, Meyerson M, Lander ES.
Nature. 2010 Jan; 6(1):99-103
}
\examples{
data(seqsumm_HCC1143)
head(seqsumm_HCC1143)
}
\keyword{datasets}
