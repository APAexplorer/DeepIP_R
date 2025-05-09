% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeepIP_funclib.R
\name{kcountByPos}
\alias{kcountByPos}
\title{kcountByPos counts kgrams at each positions in a fasta file.}
\usage{
kcountByPos(
  fafile = NULL,
  grams = NULL,
  k = NULL,
  from = NA,
  to = NA,
  sort = TRUE,
  topn = 50,
  perc = FALSE
)
}
\arguments{
\item{fafile}{a fa file name}

\item{grams}{vector of grams}

\item{k}{integer, e.g., 6}

\item{from}{if NA, then from=1}

\item{to}{if NA, then to=min length in fafile}

\item{sort}{if TRUE, then sort grams in output table}

\item{topn}{if not NULL, then return topn grams with max counts}

\item{perc}{if TRUE, to output percentage rather than counts}
}
\value{
a count or perc table with rows being grams and cols being positions from [`from` .. (`to`-`K`+1)].
}
\description{
kcountByPos counts kgrams at each positions in a fasta file.
}
\examples{
\dontrun{
cnts=kcountByPos(fafile='test.fa', grams=c('GCCACA','GAAGTC'),
                 from=90, to=110, sort=T, topn=50, perc=FALSE)
cnts=kcountByPos(fafile='test.fa', k=6, from=90, to=110,
                 sort=T, topn=50, perc=FALSE)
cnts=kcountByPos(fafile='test.fa', grams=movAPA:::getVarGrams('MM'),
                 from=90, to=110, sort=T, topn=50, perc=TRUE)
}
}
