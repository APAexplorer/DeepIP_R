% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeepIP_funclib.R
\name{sampleNseqs}
\alias{sampleNseqs}
\title{sampleNseqs randomly select N seqs from given faFiles and output each file.}
\usage{
sampleNseqs(
  faFiles,
  N,
  sufix = NULL,
  ofiles = NULL,
  outputDir = NULL,
  addLabel = NULL,
  seed = NULL
)
}
\arguments{
\item{faFiles}{fa file names.}

\item{N}{number of sequences to be sampled.}

\item{sufix}{output file name be like <faFiles><_Ns>.fa, if NULL, then sufix is like '_5000s'.}

\item{ofiles}{must the same length as faFiles. If provided, then will ignore sufix.}

\item{outputDir}{if NULL, then output each file in the dir of corresponding `faFiles`.}

\item{addLabel}{if NULL, then add <addLabel> to seq title, like <title>:1.}

\item{seed}{random seed, if NULL then is the system time.}
}
\value{
output file names.
}
\description{
sampleNseqs randomly select N seqs from given faFiles and output each file.
}
\examples{
\dontrun{
sampleNseqs(files, N=5000, outputDir='select5000sdir')
}
}
