% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cprobit.R
\name{profile_llh}
\alias{profile_llh}
\title{Inpernal function: profile log-likelihood of lambda}
\usage{
profile_llh(lambda, beta_c, y1, y2, design_mat_diff)
}
\arguments{
\item{lambda}{The Box-Cox transformation parameter. Default is \code{NA},
indicating no need for a transformation. See \code{Details}.}

\item{beta_c}{Numeric vector of Step 1 estimates.}

\item{y1}{Numeric vector of the observed outcome at observation time 1.}

\item{y2}{Numeric vector of the observed outcome at observation time 2.}

\item{design_mat_diff}{Numeric matrix of the design matrix for difference.}
}
\value{
Returns the profile log likelihood (not the negative value).
}
\description{
Inpernal function: profile log-likelihood of lambda
}
