

###############################################################################
#validMDBer_Types function definition
#	Calculates the derivative of the valid moment conditions vector for one subject.  Allows unbalanced data.
#
#Input:	
#	yvec		vector of responses.  organized by subject, by time within subject
#	subjectIndex	location of subject i within yvec
#	Zmat		design matrix for time-independent covariates ((N*T) x K0)
#	Xmat		design matrix for time-dependent covariates ((N*T) x Ktv)
#	covTypeVec	vector indicating the type of each time-dependent covariate
#	betaI		estimates of the model parameters (1+K0+Ktv x 1)
#	N		number of subjects
#	T		number of time points, subject i
#	Tmax		maximum number of time points observed

#Output:
#	dBetagEst_i	derivative of vector of valid moment conditions for i^th subject
###############################################################################

#' Generalized Method of Moments Valid Moment Combination Derivatives for one Subject of Longitudinal Bernoulli Responses, User-Defined Types
#' 
#' This function calculates the values of the derivatives of all valid moment conditions for a single subject in a longitudinal study with binary outcomes.  It allows for unbalanced data, and uses user-defined types of time-dependent covariates to determine validity of moment conditions.  The function return a matrix of derivatives for all valid moment condition for subject i.  
#' @param yvec The vector of responses, ordered by subject, time within subject.
#' @param subjectIndex The location of the first index of subject i responses within yvec.  
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate.  
#' @param betaI The current or initial estimates of the model parameters.  
#' @param N The number of subjects.  
#' @param T The number of time points for subject i.  
#' @param Tmax The maximum number of times of observation among all subjects.  
#' @keywords GMM
#' @export
#' @examples
#' validMDBer_Types()

validMDBer_Types = function(yvec,subjectIndex,Zmat,Xmat,covTypeVec,betaI,T,Tmax){

####################
# DEFINE CONSTANTS #
####################

if(!is.matrix(Zmat)){K0 = 0}
else if(is.matrix(Zmat)){K0 = ncol(Zmat)}
Ktv = ncol(Xmat)
K = 1+K0+Ktv # TOTAL NUMBER OF PARAMETERS #

K1 = 0
K2 = 0
K3 = 0
K4 = 0

for(k in 1:Ktv)
{
	if(covTypeVec[k]==1){K1 = K1+1}
	if(covTypeVec[k]==2){K2 = K2+1}
	if(covTypeVec[k]==3){K3 = K3+1}
	if(covTypeVec[k]==4){K4 = K4+1}
}

# MAXIMUM NUMBER OF ESTIMATING EQUATIONS (PER SUBJECT) #
Lmax = 1*Tmax + K0*Tmax  + (Tmax^2)*K1 + Tmax*(Tmax+1)/2*K2 + Tmax*K3 + Tmax*(Tmax+1)/2*K4

####################
####################



##################################
# CALCULATE VALUES FOR SUBJECT i #
##################################

yvec_i = yvec[subjectIndex:(subjectIndex+T)]

# CALCULATE MEAN AND SYSTEMATIC ESTIMATES FOR SUBJECT i #
mu_i = rep(0,T)
eta_i = rep(0,T)

for(t in 1:T)
{
	# PULL PREDICTORS FOR SUBJECT i, TIME t #
	if(K0!=0){zmat_it = Zmat[subjectIndex+t-1,]}
	xmat_it = Xmat[subjectIndex+t-1,]
	if(K0==0){zx_it = c(1,xmat_it)}
	else if(K0!=0){zx_it = c(1,zmat_it,xmat_it)}

	# NOTE: THIS IS SPECIFIC TO THE BERNOULLI #
	eta_i[t] = zx_it %*% betaI
	mu_i[t] = exp(eta_i[t])/(1+exp(eta_i[t]))
}

##################################
##################################







#########################################################
#	DEFINE DERIVATIVE MATRICES FOR SUBJECT i	#
#		UNIQUE TO BERNOULLI!!			##########################################################


##	FIRST DERIVATIVE OF MEAN	##

dBetamu_i = matrix(0,T,K)
for(t in 1:T)
{
dCount = 1

	# INTERCEPT #
	dBetamu_i[t,dCount] = (1)*mu_i[t]*(1-mu_i[t])
	dCount = dCount+1

	# TIC #
	if(K0!=0)
	{
	for(j in 1:K0)
	{
		dBetamu_i[t,dCount] = (Zmat[subjectIndex+t-1,j])*mu_i[t]*(1-mu_i[t])
		dCount = dCount+1
	}
	}

	# TDC #
	for(j in 1:ncol(Xmat))
	{
		dBetamu_i[t,dCount] = (Xmat[subjectIndex+t-1,j])*mu_i[t]*(1-mu_i[t])
		dCount = dCount+1
	}
}

##	PART OF SECOND DERIVATIVE OF MEAN	##

d2Betamu_i_part = matrix(0,T,K)
for (t in 1:T)
{
	for (k in 1:K)
	{
		d2Betamu_i_part[t,k] = dBetamu_i[t,k]*(1-2*mu_i[t])
	}
}




#################################################
#	MATRIX OF DERIVATIVES FOR SUBJECT i	#
#		GENERAL!!			##################################################

dBetag_i = matrix(0,Lmax,K)
count = 1

## INTERCEPT EE'S: IDENTICAL FOR ALL TIMES ##
for(t in 1:T)
{
	s = t
	j = 1
	for(k in 1:K)
	{
		dBetag_i[count,k] = (-1)*dBetamu_i[s,j]*dBetamu_i[t,k] + (1)*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
	}
	count = count+1
}
# UPDATE COUNT FOR MISSING TIMES #
count = count + (Tmax-T)


## TIC EE'S: IDENTICAL FOR ALL TIMES ##
if(K0!=0)
{
for(j in 1:K0)
{
	for(t in 1:T)
	{
		s=t
		for(k in 1:K)
		{
			dBetag_i[count,k] = (-1)*dBetamu_i[s,1+j]*dBetamu_i[t,k] + (Zmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
		}
		count = count+1
	}
	# UPDATE COUNT FOR MISSING TIMES #
	count = count + (Tmax-T)
}
}

## TDC EE'S: CONDITION ON TYPE OF TDC ##
for (j in 1:Ktv)
{
 	if(covTypeVec[j]==1) #TYPE I
	{
		for (s in 1:T)
		{
			for (t in 1:T)
			{
				for(k in 1:K)
				{
					dBetag_i[count,k] = (-1)*dBetamu_i[s,1+K0+j]*dBetamu_i[t,k] + (Xmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
				}
				count = count + 1
			}
			# UPDATE COUNT FOR MISSING TIMES #
			count = count + (Tmax-T);
		}
		# UPDATE COUNT FOR MISSING TIMES #
		count = count + Tmax*(Tmax-T);
	}
	
	else if(covTypeVec[j]==2) # TYPE II
	{
		for (s in 1:T)
		{
			for (t in 1:s)
			{
				for(k in 1:K)
				{
					dBetag_i[count,k] = (-1)*dBetamu_i[s,1+K0+j]*dBetamu_i[t,k] + (Xmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
				}
				count = count + 1
			}
		}
		# UPDATE COUNT FOR MISSING TIMES #
		count = count + (1/2)*(Tmax*(Tmax+1)-T*(T+1))
	}

	else if(covTypeVec[j]==3) # TYPE III
	{
		for (s in 1:T)
		{
			t = s
			for(k in 1:K)
			{
				dBetag_i[count,k] = (-1)*dBetamu_i[s,1+K0+j]*dBetamu_i[t,k] + (Xmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
			}
			count = count + 1
		}
		# UPDATE COUNT FOR MISSING TIMES #		
		count = count + (Tmax-T)
	}

	else # TYPE IV
	{
		for (t in 1:T)
		{
			for (s in 1:t)
			{
				for(k in 1:K)
				{
					dBetag_i[count,k] = (-1)*dBetamu_i[s,1+K0+j]*dBetamu_i[t,k] + (Xmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
				}
				count = count + 1
			}
		}
		# UPDATE COUNT FOR MISSING TIMES #
		count = count + (1/2)*(Tmax*(Tmax+1)-T*(T+1))
	}
} # END TDC EE'S #

# dBetag_i IS NOW A (L x K) MATRIX OF DERIVATIVES FOR SUBJECT i #

dBetag_i
} # END validMDBer_Types #






