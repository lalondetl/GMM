

#' Two-Step Generalized Method of Moments, Truncated Count Component of Longitudinal Hurdle Model
#'
#' This function calculates the Generalized Method of Moments (GMM) parameter estimates and standard errors for the zero-truncated count component ("positive Poisson") of a hurdle model for longitudinal excess zero count responses.  This is modeled similarly to a Positive Poisson Zero-Truncated Regression using a log link.  The function allows for unbalanced data, meaning subjects can have different numbers of times of observation.  Both time-independent covariates and time-dependent covariates can be accommodated.  Time-dependent covariates can be handled either by specifying the type of each time-dependent covariate, or by allowing the data to determine appropriate moment conditions through the extended classification method.  Data must be organized by subject, and an intercept term is assumed.  The function outputs a list with parameter estimates betaHat along with parameter covariance estimates covEst.  
#' @param y	The vector of positive count responses.  This vector must be organized by subject, and by time within subject ((sum(Tvec)) x 1).  
#' @param subjectID	The vector of subject ID values for each response ((sum(Tvec)) x 1).  
#' @param Zmat The design matrix for time-independent covariates ((sum(Tvec)) x K0).  
#' @param Xmat The design matrix for time-dependent covariates ((sum(Tvec)) x Ktv).  
#' @param Tvec The vector of times for each subject.  
#' @param N The number of subjects.  
#' @param mc The method of identifying appropriate moment conditions, either 'EC' for extended classification (default) or 'Types' for user-identified types.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate, according to the order of the columns of Xmat.  
#' @param betaI	The initial parameter estimates.  
#' r_c The vector of residuals from hurdle GEE with independent working correlation structure.  
#' @keywords GMM
#' @export
#' @examples
#' TSGMM_PP()


TSGMM_PP = function(yvec,subjectID,Zmat,Xmat,Tvec,N,mc='EC',covTypeVec=c(-1),betaI=c(1),r_c=c(-1)){


####################
# DEFINE CONSTANTS #
####################

if(!is.matrix(Zmat)){K0 = 0}
else if(is.matrix(Zmat)){K0 = ncol(Zmat)}
Ktv = ncol(Xmat)
K = 1+K0+Ktv # TOTAL NUMBER OF PARAMETERS #

Tmax = max(Tvec)

K1 = 0
K2 = 0
K3 = 0
K4 = 0

if(covTypeVec[1] != -1)
{
for(k in 1:Ktv)
{
	if(covTypeVec[k]==1){K1 = K1+1}
	if(covTypeVec[k]==2){K2 = K2+1}
	if(covTypeVec[k]==3){K3 = K3+1}
	if(covTypeVec[k]==4){K4 = K4+1}
}
}

####################
####################

# CONSTRUCT betaI AND r_c IF NECESSARY #
if(betaI == c(-1)){betaI = rep(0,(1+ncol(Zmat)+ncol(Xmat)))}




##########################################################
## MAXIMUM NUMBER OF ESTIMATING EQUATIONS (PER SUBJECT) ##
##########################################################

if(mc=='Types'){Lmax = 1*Tmax + K0*Tmax  + (Tmax^2)*K1 + Tmax*(Tmax+1)/2*K2 + Tmax*K3 + Tmax*(Tmax+1)/2*K4}

##########################################
# IF NECESSARY, FIND TYPES OF PARAMETERS #
##########################################

if(mc=='EC')
{
	alpha = 0.05/(Tmax^2)
	types = validComb_EC_PP(yvec,Zmat,Xmat,betaI,Tvec,alpha,r_c) 
	Lmax = 1*Tmax + K0*Tmax  + sum(types)
}

##########################################################
##########################################################



#############################################################
# DEFINE QUADRATIC FORM FUNCTION TAKING ONLY BETAI AS INPUT #
#############################################################

QuadForm = function(beta){
	G = rep(0,Lmax)
	VN = matrix(0,Lmax,Lmax)
	Count = rep(0,Lmax) # Vector of denominators for G and VN #

	# FILL G AND VN USING VALIDMCPP FOR EACH SUBJECT #
	for(i in 1:N)
	{
		subjectIndex = sum(Tvec[0:(i-1)])+1
		if(mc=='EC'){Est_i = validMCPP_EC(yvec,subjectIndex,Zmat,Xmat,beta,Tvec[i],Tmax,Count,types)}
		if(mc=='Types'){Est_i = validMCPP_Types(yvec,subjectIndex,Zmat,Xmat,covTypeVec,beta,Tvec[i],Tmax,Count)}

		gEst_i = Est_i[[1]]
		Count = Est_i[[2]]

		G = G + gEst_i
		VN = VN + gEst_i%*%t(gEst_i)
	}
	
	G = G / Count

	# CREATE DIVISOR MATRIX D FROM COUNT #
	D = matrix(0,Lmax,Lmax)
	for(i in 1:Lmax)
	{
		for(j in 1:Lmax)
		{
			D[i,j] = min(Count[i],Count[j])
		}
	}
	W = MASS::ginv(VN / D)

	# QUADRATIC FUNCTION TO BE OPTIMIZED #
	QF = t(G) %*% W %*% G

	QF
} # END QUADFORM #

#############################################################
#############################################################





######################################################################################
# GMM COEFFICIENTS ARE OBTAINED BY MINIMIZING QUADFORM, WITH BETAI AS INITIAL VALUES #
######################################################################################

betahat = optim(betaI, QuadForm)$par

######################################################################################
######################################################################################





########################################################################
# VARIANCE ESTIMATE IS OBTAINED USING DERIVATIVES WITH RESPECT TO BETA #
########################################################################


dBetaG = matrix(0,Lmax,K)
VN = matrix(0,Lmax,Lmax)
Count = rep(0,Lmax) # Vector of denominators for G and VN #

for(i in 1:N)
{
	subjectIndex = sum(Tvec[0:(i-1)])+1

	if(mc=='EC'){Est_i = validMCPP_EC(yvec,subjectIndex,Zmat,Xmat,betahat,Tvec[i],Tmax,Count,types)}
	if(mc=='Types'){Est_i = validMCPP_Types(yvec,subjectIndex,Zmat,Xmat,covTypeVec,betahat,Tvec[i],Tmax,Count)}

	gEst_i = Est_i[[1]]
	Count = Est_i[[2]]
	VN = VN + gEst_i%*%t(gEst_i)

	if(mc=='EC'){dBetagEst_i = validMDPP_EC(yvec,subjectIndex,Zmat,Xmat,betahat,Tvec[i],Tmax,types)}
	if(mc=='Types'){dBetagEst_i = validMDPP_Types(yvec,subjectIndex,Zmat,Xmat,covTypeVec,betahat,Tvec[i],Tmax)}

	dBetaG = dBetaG + dBetagEst_i
}

# CREATE DIVISOR MATRIX D FROM COUNT #
D = matrix(0,Lmax,Lmax)
for(i in 1:Lmax)
{
	for(j in 1:Lmax)
	{
		D[i,j] = min(Count[i],Count[j])
	}
}

## CREATE DIVISOR MATRIX FOR DBETAG FROM COUNT ##
Divisor = matrix(c(rep(Count,K)),length(Count),K)

dBetaG = dBetaG / Divisor
W = MASS::ginv(VN / D)

AsymptoticWeight = t(dBetaG) %*% W %*% dBetaG
AsymptoticCovariance = (1/N)*MASS::ginv(AsymptoticWeight)

########################################################################
########################################################################


list(betaHat=betahat, covEst = AsymptoticCovariance)

} # END TSGMM_PP #


