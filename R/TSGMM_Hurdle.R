

#' Two-Step Generalized Method of Moments, Longitudinal Hurdle Model
#'
#' This function calculates the Generalized Method of Moments (GMM) parameter estimates and standard errors for longitudinal excess zero count responses using a hurdle model formulation.  The function allows for unbalanced data, meaning subjects can have different numbers of times of observation.  Both time-independent covariates and time-dependent covariates can be accommodated.  Time-dependent covariates can be handled either by specifying the type of each time-dependent covariate, or by allowing the data to determine appropriate moment conditions through the extended classification method.  Data must be organized by subject, and an intercept term is assumed for both the logistic and truncated count components of the model.  The function outputs a list with parameter estimates betaHat_l and betaHat_c for the logistic and count components, respectively, along with parameter covariance estimates covEst_l and covEst_c.  
#' @param y	The vector of response counts, ordered by subject, time within subject.  
#' @param subjectID	The vector of subject ID values for each response.  
#' @param  N The number of subjects.  
#' @param Tvec The vector of times for each subject (N x 1) 
#' @param X_l The design matrix for time-dependent covariates in the logistic component of the model, excluding an intercept.  
#' @param Z_l The design matrix for time-independent covariates in the logistic component of the model, excluding an intercept. 
#' @param X_c The design matrix for time-dependent covariates in the truncated count component of the model, excluding an intercept. 
#' @param Z_c The design matrix for time-independent covariates in the truncated count component of the model, excluding an intercept. 
#' @param mc The method of identifying appropriate moment conditions, either 'EC' for extended classification (default) or 'Types' for user-identified types.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate, according to the order of the columns of X_l and X_c.  
#' @keywords GMM
#' @export
#' @examples
#' TSGMM_Hurdle()



TSGMM_Hurdle = function(y,subjectID,N,Tvec,X_l,Z_l,X_c,Z_c,mc='EC',covTypeVec_l=c(-1),covTypeVec_c=c(-1)){


# CREATE ONE DESIGN FOR GEE #
if(!is.matrix(Z_l)){ZX_l = cbind(rep(1,nrow(X_l)),X_l)}
if(is.matrix(Z_l)){ZX_l = cbind(rep(1,nrow(X_c)),Z_l,X_l)}

if(!is.matrix(Z_c)){ZX_c = cbind(rep(1,nrow(X_c)),X_c)}
if(is.matrix(Z_c)){ZX_c = cbind(rep(1,nrow(X_c)),Z_c,X_c)}

# APPLY  HURDLE GEE #
HGEE_output = HurdleGEE(y,subjectID,N,ZX_l,ZX_c)

betaI_l = HGEE_output$betaHat_l
betaI_c = HGEE_output$betaHat_c

r_l = HGEE_output$r_l
r_c = HGEE_output$r_c


# CALL ADJUSTED VERSION OF TSGMM_BER #
TSGMM_Hurdle_l = TSGMM_c0(y,subjectID,Z_l,X_l,Tvec,N,mc,covTypeVec_l,betaI_l,r_l)

# CALL ADJUSTED VERSION OF TSGMM_PP #
TSGMM_Hurdle_c = TSGMM_PP(y,subjectID,Z_c,X_c,Tvec,N,mc,covTypeVec_c,betaI_c,r_c)

# COMBINE #
betaHat_l = TSGMM_Hurdle_l$betaHat
betaHat_c = TSGMM_Hurdle_c$betaHat

covEst_l = TSGMM_Hurdle_l$covEst
covEst_c = TSGMM_Hurdle_c$covEst

list(betaHat_l=betaHat_l, betaHat_c=betaHat_c, covEst_l=covEst_l, covEst_c=covEst_c)

} # END PARENT FUNCTION #



