

#' Two-Step Generalized Method of Moments, Longitudinal Outcome
#'
#' This function calculates the Generalized Method of Moments (GMM) parameter estimates and standard errors for longitudinal data.  The function allows for unbalanced data, meaning subjects can have different numbers of times of observation.  Both time-independent covariates and time-dependent covariates can be accommodated.  Time-dependent covariates can be handled either by specifying the type of each time-dependent covariate, or by allowing the data to determine appropriate moment conditions through the extended classification method.  
#' @param yvec The vector of responses, ordered by subject, time within subject.
#' @param subjectID The vector of subject ID values for each response.  
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param Tvec The vector of times for each subject.  
#' @param N The number of subjects.  
#' @param mc The method of identifying appropriate moment conditions, either 'EC' for extended classification (default) or 'Types' for user-identified types.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate, according to the order of the columns of Xmat.  
#' @param family The family indicates the response type, and follows standard R family notation.  
#' @keywords GMM
#' @export
#' @examples
#' TSGMM()




TSGMM = function(yvec,subjectID,Zmat,Xmat,Tvec,N,mc='EC',covTypeVec=c(-1),family='gaussian'){

# CONDITION ON FAMILY, CALL APPROPRIATE TSGMM FUNCTION #

if(family=='gaussian'){Model=TSGMM_Nor(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(family=='binomial'){Model=TSGMM_Ber(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(family=='poisson'){Model=TSGMM_Poi(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(family=='ppoisson'){Model=TSGMM_PP(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

#if(family=='beta'){Model=TSGMM_Beta(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}


# HURDLE AND ZIP NEED TO BE HANDLED SEPARATELY BECAUSE OF THE JOINT MODELS #

#if(family=='hurdle'){Model=TSGMM_Hurdle(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

#if(family=='zip'){Model=TSGMM_ZIP(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}


list(betaHat=Model$betahat, covEst = Model$covEst)

} # END TSGMM #



