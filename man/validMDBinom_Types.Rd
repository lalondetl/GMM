% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validMDBinom_Types.R
\name{validMDBinom_Types}
\alias{validMDBinom_Types}
\title{Generalized Method of Moments Valid Moment Combination Derivatives for one Subject of Longitudinal Count (number of events from n trials) Responses, User-Defined Types}
\usage{
validMDBinom_Types(ymat, subjectIndex, Zmat, Xmat, covTypeVec, betaI, T, Tmax)
}
\arguments{
\item{ymat}{The matrix of responses, ordered by subject, time within subject.  The first column is the number of successes, the second the number of failures.}

\item{subjectIndex}{The location of the first index of subject i responses within ymat.}

\item{Zmat}{The design matrix for time-independent covariates.}

\item{Xmat}{The design matrix for time-dependent covariates.}

\item{covTypeVec}{The vector indicating the type of each time-dependent covariate.}

\item{betaI}{The current or initial estimates of the model parameters.}

\item{T}{The number of time points for subject i.}

\item{Tmax}{The maximum number of times of observation among all subjects.}

\item{N}{The number of subjects.}
}
\description{
This function calculates the values of the derivatives of all valid moment conditions for a single subject in a longitudinal study with count (0-n) outcomes.  It is assumed that the count represents the number of events from n identical trials, and that n is equal for all subjects and times.  This is modeled similarly to a Logistic Regression for Binomial responses.  It allows for unbalanced data, and uses user-defined types of time-dependent covariates to determine validity of moment conditions.  The function returns a matrix of derivatives for all valid moment condition for subject i.
}
\examples{
validMDBinom_Types()
}
\keyword{GMM}
