% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeepIP_funclib.R
\name{getTxdbIntrons}
\alias{getTxdbIntrons}
\title{getTxdbIntrons gets intron ranges in TXDB.}
\usage{
getTxdbIntrons(txdb, chrs = NULL)
}
\arguments{
\item{txdb}{a TXDB object , or a data.frame with at least columns: seqnames, start, end ,strand, to get intron ranges.}

\item{chrs}{a string vector storing chrs to subset, default is NULL.
If provided, will try to format the seqnames in txdb to make consistent with chrs.}
}
\value{
a GRanges of introns
}
\description{
getTxdbIntrons gets intron ranges in TXDB.
}
