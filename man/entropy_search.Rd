% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Entropy search.R
\name{entropy_search}
\alias{entropy_search}
\title{Entropy Search}
\usage{
entropy_search(
  X,
  f,
  N,
  lower = 0,
  upper = 1,
  M = 100,
  K = 1000,
  particles = 5000,
  g = NULL,
  update.particles = 5,
  temps.start = NULL,
  walk.sd = NULL,
  discrete.type = "EI",
  prop.EI = 0.8,
  direction = "both",
  updateGP = 1,
  track.progress = TRUE,
  rePT = TRUE
)
}
\arguments{
\item{walk.sd}{}
}
\description{
Entropy Search
}
