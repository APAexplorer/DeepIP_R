% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DeepIP_funclib.R
\name{statDeepResByGrams}
\alias{statDeepResByGrams}
\title{statDeepResByGrams gets metrics for each gram for a .csv prediction file from DeepIP_test.py.}
\usage{
statDeepResByGrams(deepResCsv, ofile = NULL, verbose = FALSE, plot = TRUE)
}
\arguments{
\item{deepResCsv}{one csv files from DeepIP_test.py, with the title column like ">PA615:chr1;+;16972964;AACAAA:1".
Or if the `gram` column is already in the csv, then will be used instead of parsing the title column.}

\item{ofile}{output .csv or other file to store the metrics table. If it is NULL, then not output.}

\item{verbose}{TRUE to print the metrics table.}

\item{plot}{TRUE to plot F1 ROC Sn Sp to compare grams.}
}
\value{
A dataframe of the metrics, with each column one gram in file of `deepResCsvs` and each row a metric.
Metrics are: TP_11	FP_01	TN_00	FN_10	Sensitivity	Specificity	Precision	Recall	F1	ROC	AUC.
grams are obtained from the last chars of the seq title like ">PA615:chr1;+;16972964;AACAAA:1".
In the bar plot, the grams are ordered by F1 desc. AATAAA, ATTAAA, and NOPAS are highlighted in darkgrey.
}
\description{
statDeepResByGrams gets metrics for each gram for a .csv prediction file from DeepIP_test.py.
}
\examples{
\dontrun{
stats=statDeepResByGrams('train.10000T.10000F.1.epoch100_ON_test.all.csv',
                         ofile=NULL,
                         plot=TRUE)
}
}
