

#' Generalized Method of Moments Valid Moment Combinations for one Subject of Longitudinal Count (number of events from n trials) Responses, using Extended Classification
#' 
#' This function calculates the values of the valid moment conditions for a single subject in a longitudinal study with count (0-n) outcomes.  It is assumed that the count represents the number of events from n identical trials, and that n is equal for all subjects and times.  This is modeled similarly to a Logistic Regression for Binomial responses.  It allows for unbalanced data, and uses the extended classification method to determine validity of moment conditions.  The function returns a vector of valid moment condition values for subject i, along with an updated count vector of the number of valid moment conditions.  
#' @param ymat The matrix of responses, ordered by subject, time within subject.  The first column is the number of successes, the second the number of failures.  
#' @param subjectIndex The location of the first index of subject i responses within ymat.  
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param betaI The current or initial estimates of the model parameters.  
#' @param T The number of time points for subject i.  
#' @param Tmax The maximum number of times of observation among all subjects.  
#' @param Count A vector of running counts of the number of valid moment conditions for all subjects.  
#' @param types A matrix indicating valid moment conditions for all time-dependent covariate parameters.  
#' @keywords GMM
#' @export
#' @examples
#' validMCBinom_EC()

validMCBinom_EC = function(ymat,subjectIndex,Zmat,Xmat,betaI,T,Tmax,Count,types){

####################
# DEFINE CONSTANTS #
####################

if(!is.matrix(Zmat)){K0 = 0}
else if(is.matrix(Zmat)){K0 = ncol(Zmat)}
Ktv = ncol(Xmat)
K = 1+K0+Ktv # TOTAL NUMBER OF PARAMETERS #

Lmax = 1*Tmax + K0*Tmax + sum(types)

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

#Time-dependent covariate moments: condition on EC
for (k in 1:Ktv)
{
	types_k = types[,((Tmax*(k-1)+1):(Tmax*(k-1)+Tmax))]

	for(s in 1:T){
		for(t in 1:T){
			if(types_k[s,t] == 1){
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(ymat_i[t]-mu_i[t])
				Count[count] = Count[count]+1
				count = count+1
			}
		}
	
		# UPDATE COUNT FOR MISSING TIMES #
		if(T<Tmax){count = count + sum(types_k[s,((T+1):Tmax)])}
	}

	# UPDATE COUNT FOR MISSING TIMES #
	if(T<Tmax){count = count + sum(types_k[((T+1):Tmax),])}

} #end TIME-DEPENDENT COVARIATE MOMENTS

#gEst_i is now a (L x 1) vector of valid moment conditions for subject i
#Count has been incremented for each valid moment condition of subject i

list(gEst_i,Count)

} # end validMCBinom_EC #







