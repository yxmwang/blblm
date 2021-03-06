% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblog.R
\name{predict.blblog}
\alias{predict.blblog}
\title{the prediction interval in blblog}
\usage{
\method{predict}{blblog}(object, new_data, confidence = FALSE, level = 0.95, ...)
}
\arguments{
\item{object}{fitted model from blblog}

\item{new_data}{a new data frame}

\item{confidence}{boolean}

\item{level}{numeric, the confidence interval level}

\item{...}{additional arguements/inputs}
}
\value{
the confidence interval of the estimates
}
\description{
the prediction interval in blblog
}
