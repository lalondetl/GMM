% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validMCBinom_Types.R
\name{validMCBinom_Types}
\alias{validMCBinom_Types}
\title{Generalized Method of Moments Valid Moment Combinations for Longitudinal Count (number of events from n trials) Responses, User-Defined Types}
\usage{
validMCBinom_Types(ymat, subjectIndex, Zmat, Xmat, covTypeVec, betaI, T, Tmax,
  Count)
}
\arguments{
\item{ymat}{The matrix of responses, ordered by subject, time within subject.  The first column is the number of successes, the second the number of failures.}

\item{Zmat}{The design matrix for time-independent covariates.}

\item{Xmat}{The design matrix for time-dependent covariates.}

\item{covTypeVec}{The vector indicating the type of each time-dependent covariate.}

\item{betaI}{The current or initial estimates of the model parameters.}

\item{T}{The number of time points for subject i.}

\item{Tmax}{The maximum number of times of observation among all subjects.}

\item{Count}{A vector of running counts of the number of valid moment conditions for all subjects.}
}
\description{
This function calculates the values of valid moment combinations for two-step Generalized Method of Moments with user-defined types of time-dependent covariates, applied to longitudinal data with count (0-n) outcomes.  It is assumed that the count represents the number of events from n identical trials, and that n is equal for all subjects and times.  This is modeled similarly to a Logistic Regression for Binomial responses.  It allows for unbalanced longitudinal data, meaning subjects can be observed for different numbers of times.  The function returns a vector "types" indicating validity of different moments conditions.
}
\examples{
validMCBinom_Types()
}
\keyword{GMM}
