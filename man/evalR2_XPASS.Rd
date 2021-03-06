% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_XPASS.R
\name{evalR2_XPASS}
\alias{evalR2_XPASS}
\title{Evaluate the approximated predictive R2 using an independent GWAS summary statistics.}
\usage{
evalR2_XPASS(pm, file_z_pred, file_predgeno)
}
\arguments{
\item{pm}{The output `mu` obtained from the `XPASS` function, which stores the posterior means computed by XPASS.}

\item{file_z_pred}{Summary statistics file of the test data. The format is the same as XPASS input.}

\item{file_predgeno}{Prefix of the reference genotypes in plink format.}
}
\value{
Approximated predictive R^2.
}
\description{
This function approximates the predictive R2 of XPASS using an independent GWAS summary statistics.
}
\examples{
See vignette.
}
