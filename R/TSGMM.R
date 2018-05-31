

#' Two-Step Generalized Method of Moments, Longitudinal Outcome
#'
#' This function calculates the Generalized Method of Moments (GMM) parameter estimates and standard errors for longitudinal data.  The function allows for unbalanced data, meaning subjects can have different numbers of times of observation.  Both time-independent covariates and time-dependent covariates can be accommodated.  Time-dependent covariates can be handled either by specifying the type of each time-dependent covariate, or by allowing the data to determine appropriate moment conditions through the extended classification method.  While the GMM estimation method does not require a data-generating process, and therefore response are not assumed to be distributed according to known distributions, data response types can loosely be associated with traditional response distributions.  These include continuous outcomes (similar to Normal models), binary outcomes (similar to Bernoulli logistic models), proportions (similar to Beta regression), count events (similar to Binomial logistic models), counts (similar to Poisson regression), positive-counts (similar to Positive Poisson regression), hurdle outcomes (similar to Hurdle Poisson regression), and zero-inflated outomes (similar to Zero-Inflated Poisson models).  
#' @param yvec The vector of responses, ordered by subject, time within subject.  If the response type is 'count-events' (similar to Binomial data) the response is assumed to be a matrix with two columns: the first column will be the count of successes (or events), and the second column will be the count of failures (or non-events).  
#' @param subjectID The vector of subject ID values for each response.  
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param Tvec The vector of times for each subject.  
#' @param N The number of subjects.  
#' @param mc The method of identifying appropriate moment conditions, either 'EC' for extended classification (default) or 'Types' for user-identified types.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate, according to the order of the columns of Xmat.  
#' @param response The response type, either 'continuous', 'binary', 'proportion', 'count-events', 'count', 'positive-count', 'hurdle', or 'zi'.  
#' @keywords GMM
#' @export
#' @examples
#' TSGMM()




TSGMM = function(yvec,subjectID,Zmat,Xmat,Tvec,N,mc='EC',covTypeVec=c(-1),response='continuous'){

# CONDITION ON RESPONSE, CALL APPROPRIATE TSGMM FUNCTION #

if(type=='continuous'){Model=TSGMM_Nor(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(type=='binary'){Model=TSGMM_Ber(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(type=='proportion'){Model=TSGMM_Beta(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(type=='count-events'){Model=TSGMM_Binom(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(type=='count'){Model=TSGMM_Poi(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

if(type=='positive-count'){Model=TSGMM_PP(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}



# HURDLE AND ZIP NEED TO BE HANDLED SEPARATELY BECAUSE OF THE JOINT MODELS #

#if(type=='hurdle'){Model=TSGMM_Hurdle(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}

#if(type=='zi'){Model=TSGMM_ZI(yvec,subjectID,Zmat,Xmat,Tvec,N,mc,covTypeVec)}


list(betaHat=Model$betahat, covEst = Model$covEst)

} # END TSGMM #



