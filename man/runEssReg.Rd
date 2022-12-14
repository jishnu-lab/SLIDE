% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runEssReg.R
\name{runEssReg}
\alias{runEssReg}
\title{Run Essential Regression for cross-validation.}
\usage{
runEssReg(
  train_y,
  train_x_raw,
  train_x_std,
  valid_x,
  delta,
  lambda,
  thresh_fdr,
  rep_cv,
  alpha_level
)
}
\arguments{
\item{train_y}{a vector of numeric values; the response (training set)}

\item{train_x_raw}{a matrix or data frame of numeric values; the data (training set)}

\item{train_x_std}{a matrix or data frame of numeric values; the data (training set)}

\item{valid_x}{a matrix or data frame of numeric values; the data (validation set)}

\item{delta}{\eqn{\delta}, a numerical value used for thresholding}

\item{lambda}{\eqn{\lambda}, a numerical constant used in thresholding}

\item{thresh_fdr}{a numerical constant used for thresholding the correlation matrix to
control the false discovery rate, default is 0.2}

\item{rep_cv}{number of replicates for inner cross-validation of \code{delta}}

\item{alpha_level}{\eqn{\alpha}, a numerical constant used in confidence interval calculation}
}
\value{
a list containing the training model and the predicted values of the response
using \code{valid_x} and using \code{train_x_std}
}
\description{
Runs Essential Regression to create a model and predicts response values
for the validation set.
}
