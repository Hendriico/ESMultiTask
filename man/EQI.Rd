% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/criteria.R
\name{EQI}
\alias{EQI}
\title{Expected Quantile Improvement}
\usage{
EQI(
  x,
  model,
  model.from = c("laGP", "km"),
  comp.model.noise.var = 0,
  beta = 0.9,
  qmin = NULL,
  x.design = NULL,
  predf = laGP::predGPsep
)
}
\arguments{
\item{predf}{}
}
\description{
Expected Quantile Improvement
}
