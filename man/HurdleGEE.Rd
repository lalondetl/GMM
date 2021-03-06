% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HurdleGEE.R
\name{HurdleGEE}
\alias{HurdleGEE}
\title{Generalized Estimating Equations, Longitudinal Hurdle Model}
\usage{
HurdleGEE(y, subjectID, N, X_l, X_c)
}
\arguments{
\item{y}{The vector of response counts, ordered by subject, time within subject.}

\item{subjectID}{The vector of subject ID values for each response.}

\item{N}{The number of subjects.}

\item{X_l}{The design matrix for all covariates in the logistic component of the model, including an intercept.}

\item{X_c}{The design matrix for all covariates in the truncated count component of the model, including an intercept.}
}
\description{
This function calculates the Generalized Estimating Equations (GEE) parameter estimates and standard errors for longitudinal excess zero count responses using a hurdle model formulation with independent working correlation structure, based on Dobbie and Welsh (2001).  Data must be organized by subject, and an intercept term is assumed for both the logistic and truncated count components of the model.  The function outputs a list with parameter estimates betaHat_l and betaHat_c for the logistic and count components, respectively, along with parameter covariance estimates covEst_l and covEst_c and residuals r_l and r_c.
}
\examples{
HurdleGEE()
}
\keyword{GEE}
