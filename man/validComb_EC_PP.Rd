% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validComb_EC_PP.R
\name{validComb_EC_PP}
\alias{validComb_EC_PP}
\title{Generalized Method of Moments Valid Moment Combinations for Truncated Count Component of Longitudinal Excess-Zero Responses, using Extended Classification}
\usage{
validComb_EC_PP(yvec, Zmat, Xmat, betaI, Tvec, alpha, r)
}
\arguments{
\item{yvec}{The vector of responses, ordered by subject, time within subject.}

\item{Zmat}{The design matrix for time-independent covariates ((N*T) x K0).}

\item{Xmat}{The design matrix for time-dependent covariates ((N*T) x Ktv).}

\item{betaI}{The current or initial estimates of the model parameters (1+K0+Ktv x 1).}

\item{Tvec}{The vector of times for each subject.}

\item{alpha}{The significance level for inclusion of moment conditions.}

\item{r}{The residuals from fitting an initial hurdle GEE model using the independent working correlation structure.}
}
\description{
This function calculates the values of valid moment combinations for two-step Generalized Method of Moments using the extended classification method, applied to the zero-truncated count component ("positive Poisson") of a hurdle model for excess-zero longitudinal data.  It allows for unbalanced longitudinal data, meaning subjects can be observed for different numbers of times.  The function returns a vector "types" indicating validity of different moments conditions.
}
\examples{
validComb_EC_PP()
}
\keyword{GMM}
