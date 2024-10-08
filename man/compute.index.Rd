% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_1_Index_computation_20240909.R
\name{compute.index}
\alias{compute.index}
\title{Compute the diversity index from pairwise evolutionary time differences between all nodes and tips}
\usage{
compute.index(
  time_distance_mat,
  timed_tree,
  time_window,
  metadata,
  mutation_rate,
  timescale,
  genome_length
)
}
\arguments{
\item{time_distance_mat}{Hamming distance matrix}

\item{timed_tree}{Timed tree}

\item{time_window}{Time window for the index computation}

\item{metadata}{Metadata for all tips and nodes (see details)}

\item{mutation_rate}{Mutation rate (in mutation/site/year)}

\item{timescale}{Index timescale (in years)}

\item{genome_length}{Genome length (in bp)}
}
\value{
The index of each sequence and node in the tree
}
\description{
Compute the diversity index from pairwise evolutionary time differences between all nodes and tips
}
\references{
Wirth, T, Wong, V, Vandenesch, F, Rasigade, J-P. Applied phyloepidemiology: Detecting drivers of pathogen transmission from genomic signatures using density measures. Evol Appl. 2020; 13: 1513– 1525. https://doi.org/10.1111/eva.12991
}
