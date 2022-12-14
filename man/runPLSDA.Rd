% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runPLSDA.R
\name{runPLSDA}
\alias{runPLSDA}
\title{Run partial least-squares discriminant analysis for cross-validation.}
\usage{
runPLSDA(train_y, train_x, valid_x)
}
\arguments{
\item{train_y}{a vector of numeric values; the response (training set)}

\item{train_x}{a matrix or data frame of numeric values; the data (training set)}

\item{valid_x}{a matrix or data frame of numeric values; the data (validation set)}
}
\value{
a list containing the training model and the predicted class probabilities of the response
using \code{valid_x} and \code{train_x}
}
\description{
Creates a partial least-squares discriminant analysis classifier and
predicts response values from the validation data.
}
