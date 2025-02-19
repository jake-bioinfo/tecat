% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f06-Telomere_length.R
\name{determine_threshold}
\alias{determine_threshold}
\title{determine_threshold}
\usage{
determine_threshold(
  telomere_list = NULL,
  platform = "pb",
  sample_ratio = 0.15,
  threads = 0
)
}
\arguments{
\item{telomere_list}{A list of data frames containing the telomere
sequences to be processed.}

\item{platform}{A character string specifying the sequencing platform
used to generate the data. Default is "pb" for PacBio.}

\item{sample_ratio}{A numeric value specifying the proportion of the
data to sample for threshold determination. Default is 0.15.}

\item{threads}{An integer specifying the number of parallel processors
to use. Default is 0.}
}
\value{
A list of data frames containing the start and end thresholds
for the telomere length calculation.
}
\description{
Determine the start and stop thresholds for the telomere
length calculation.
}
