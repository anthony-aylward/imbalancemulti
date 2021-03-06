% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_imbalance_data.R
\name{read_imbalance_data}
\alias{read_imbalance_data}
\title{Read Imbalance Data}
\usage{
read_imbalance_data(..., experiment_names = NULL, cores = 1)
}
\arguments{
\item{...}{paths to imbalance data files}

\item{experiment_names}{character, names of experiments}

\item{cores}{integer, number of cores to use}
}
\value{
list containing matrices of allele counts, indexed by variant id
}
\description{
Read imbalance data from multiple experiments
}
