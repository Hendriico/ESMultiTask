% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Knowledge Gradient.R
\name{KG}
\alias{KG}
\title{Knowledge Gradient}
\usage{
KG(
  Xinit,
  f,
  N = 100,
  noisy = TRUE,
  sigma = 0,
  model.type = c("laGP", "km"),
  KG.type = "Correlated",
  type = "UK",
  lower = 0,
  upper = 1,
  B = 2,
  optim.multi.start = 5
)
}
\arguments{
\item{optim.multi.start}{}
}
\description{
This algorithm uses either Correlated KG or Approximate KG criteria with Gaussian processes
to find the minimum of a blackbox function
}
