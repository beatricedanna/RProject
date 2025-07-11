% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NormalizeCPM.R
\name{NormalizeCPM}
\alias{NormalizeCPM}
\title{NormalizeCPM}
\usage{
NormalizeCPM(data_table, chunk_size = 1000)
}
\arguments{
\item{data_table}{A data.table containing counts data. It must include a column named "feature_id".}

\item{chunk_size}{Integer, optional (default = 1000). Number of columns to process at a time.}
}
\value{
A data.table with the same dimensions as data_table, where counts are normalized to CPM and log2-transformed. The "feature_id" column is mantained as the first column.
}
\description{
NormalizeCPM
}
