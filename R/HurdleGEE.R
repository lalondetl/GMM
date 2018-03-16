

#' Generalized Estimating Equations, Longitudinal Hurdle Model
#'
#' This function calculates the Generalized Estimating Equations (GEE) parameter estimates and standard errors for longitudinal excess zero count responses using a hurdle model formulation with independent working correlation structure, based on Dobbie and Welsh (2001).  Data must be organized by subject, and an intercept term is assumed for both the logistic and truncated count components of the model.  The function outputs a list with parameter estimates betaHat_l and betaHat_c for the logistic and count components, respectively, along with parameter covariance estimates covEst_l and covEst_c and residuals r_l and r_c.  
#' @param y	The vector of response counts, ordered by subject, time within subject.  
#' @param subjectID	The vector of subject ID values for each response.  
#' @param  N The number of subjects.  
#' @param X_l The design matrix for all covariates in the logistic component of the model, including an intercept.  
#' @param X_c The design matrix for all covariates in the truncated count component of the model, including an intercept. 
#' @keywords GEE
#' @export
#' @examples
#' HurdleGEE()



HurdleGEE = function(y,subjectID,N,X_l,X_c){

logistic_output = HurdleGEE_l(y,subjectID,N,X_l)
count_output = HurdleGEE_c(y,subjectID,N,X_c)

betaHat_l = logistic_output$betaHat
r_l = logistic_output$r_l
covEst_l = logistic_output$covEst

betaHat_c = count_output$betaHat
r_c = count_output$r_c 	## NOTE MISSING VALUES WHERE Y==0 ##
covEst_c = count_output$covEst

list(betaHat_l=betaHat_l,r_l=r_l,covEst_l=covEst_l,betaHat_c=betaHat_c,r_c=r_c,covEst_c=covEst_c)

} # END HurdleGEE #





