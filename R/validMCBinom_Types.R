

#' Generalized Method of Moments Valid Moment Combinations for Longitudinal Count (number of events from n trials) Responses, User-Defined Types
#' 
#' This function calculates the values of valid moment combinations for two-step Generalized Method of Moments with user-defined types of time-dependent covariates, applied to longitudinal data with count (0-n) outcomes.  It is assumed that the count represents the number of events from n identical trials, and that n is equal for all subjects and times.  This is modeled similarly to a Logistic Regression for Binomial responses.  It allows for unbalanced longitudinal data, meaning subjects can be observed for different numbers of times.  The function returns a vector "types" indicating validity of different moments conditions.  
#' @param ymat The matrix of responses, ordered by subject, time within subject.  The first column is the number of successes, the second the number of failures.  
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate.  
#' @param betaI The current or initial estimates of the model parameters.  
#' @param T The number of time points for subject i.  
#' @param Tmax The maximum number of times of observation among all subjects.  
#' @param Count A vector of running counts of the number of valid moment conditions for all subjects.  
#' @keywords GMM
#' @export
#' @examples
#' validMCBinom_Types()


validMCBinom_Types = function(ymat,subjectIndex,Zmat,Xmat,covTypeVec,betaI,T,Tmax,Count){

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

ymat_i = ymat[subjectIndex:(subjectIndex+T)]

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




##########################################################
# DEFINE VECTOR OF VALID MOMENT CONDITIONS FOR SUBJECT i #
##########################################################

gEst_i = rep(0,Lmax)
count = 1

#NOTE: gEst_i function is specific to the Bernoulli / Binomial
#Intercept term moments: identical for all times
for(t in 1:T)
{
	gEst_i[count] = (mu_i[t]/(1+exp(eta_i[t])))*(ymat_i[t]-mu_i[t])
	Count[count] = Count[count]+1
	count = count+1
}

# UPDATE COUNT FOR MISSING TIMES #
count = count + (Tmax-T)

#Time-independent covariate moments: identical for all times
if(K0!=0)
{
for(k in 1:K0)
{
	for(t in 1:T)
	{
		gEst_i[count] = (mu_i[t]/(1+exp(eta_i[t])))*Zmat[subjectIndex+t-1,k]*(ymat_i[t]-mu_i[t])
		Count[count] = Count[count]+1
		count = count+1
	}

	# UPDATE COUNT FOR MISSING TIMES #
	count = count + (Tmax-T)

}
}

#Time-dependent covariate moments: condition on type of TVC
for (k in 1:Ktv)
{
 	if(covTypeVec[k]==1) # TYPE I #
	{
		for (s in 1:T)
		{
			for (t in 1:T)
			{
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(ymat_i[t]-mu_i[t])
				Count[count] = Count[count]+1
				count = count + 1
			}

			# UPDATE COUNT FOR MISSING TIMES #
			count = count + (Tmax-T)
		}

		# UPDATE COUNT FOR MISSING TIMES #
		count = count + Tmax*(Tmax-T)
	}
	
	else if(covTypeVec[k]==2) # TYPE II #
	{
		for (s in 1:T)
		{
			for (t in 1:s)
			{
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(ymat_i[t]-mu_i[t])
				Count[count] = Count[count]+1
				count = count + 1
			}
		}

		# UPDATE COUNT FOR MISSING TIMES #
		count = count + (1/2)*(Tmax*(Tmax+1)-T*(T+1))
	}

	else if(covTypeVec[k]==3) # TYPE III #
	{
		for (s in 1:T)
		{
			gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(ymat_i[s]-mu_i[s])
			Count[count] = Count[count]+1
			count = count + 1
		}

		# UPDATE COUNT FOR MISSING TIMES #
		count = count + (Tmax-T)
	}

	else # TYPE IV #
	{
		for (t in 1:T)
		{
			for (s in 1:t)
			{
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(ymat_i[t]-mu_i[t])
				Count[count] = Count[count]+1
				count = count + 1
			}
		}

		# UPDATE COUNT FOR MISSING TIMES #
		count = count + (1/2)*(Tmax*(Tmax+1)-T*(T+1))
	}
} #end TIME-DEPENDENT COVARIATE MOMENTS

#gEst_i is now a (L x 1) vector of valid moment conditions for subject i
#Count has been incremented for each valid moment condition of subject i

list(gEst_i,Count)

} # end validMCBinom_Types #





