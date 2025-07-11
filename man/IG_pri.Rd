% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priors.R
\name{IG_pri}
\alias{IG_pri}
\title{Function to set an Inverse Gamma prior}
\usage{
IG_pri(a, b)
}
\arguments{
\item{a}{shape parameter.}

\item{b}{scale parameter.}
}
\value{
A list with the prior name and all the quantities required to implement it within the sampler.
}
\description{
Allow the user to specify an Inverse Gamma prior for a variance parameter
}
