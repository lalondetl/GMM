

#' Generalized Method of Moments Valid Moment Combinations for one Subject of Longitudinal Zero-Truncated Count Responses, using Extended Classification
#' 
#' This function calculates the values of the valid moment conditions for a single subject in a longitudinal study with zero-truncated count ("positive Poisson") outcomes.  It allows for unbalanced data, and uses the extended classification method to determine validity of moment conditions.  The function returns a vector of valid moment condition values for subject i, along with an updated count vector of the number of valid moment conditions.  
#' @param yvec The vector of responses, ordered by subject, time within subject.
#' @param subjectIndex The location of the first index of subject i responses within yvec.  
#' @param Zmat The design matrix for time-independent covariates ((N*T) x K0).  
#' @param Xmat The design matrix for time-dependent covariates ((N*T) x Ktv).  
#' @param betaI The current or initial estimates of the model parameters (1+K0+Ktv x 1).  
#' @param T The number of time points for subject i.  
#' @param Tmax The maximum number of times of observation among all subjects.  
#' @param Count A vector of running counts of the number of valid moment conditions for all subjects.  
#' @param types A matrix indicating valid moment conditions for all time-dependent covariate parameters.  
#' @keywords GMM
#' @export
#' @examples
#' validMCPP_EC()



validMCPP_EC = function(yvec,subjectIndex,Zmat,Xmat,betaI,T,Tmax,Count,types){

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

yvec_i = yvec[subjectIndex:(subjectIndex+T)]
y_ind_i = ifelse(yvec_i>0,1,0)


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

	# NOTE: THIS IS SPECIFIC TO THE POSITIVE POISSON #
	eta_i[t] = zx_it %*% betaI
	mu_i[t] = exp(eta_i[t])/(1-exp(-exp(eta_i[t])))
}

##################################
##################################






##########################################################
# DEFINE VECTOR OF VALID MOMENT CONDITIONS FOR SUBJECT i #
##########################################################


gEst_i = rep(0,Lmax)
count = 1

#NOTE: gEst_i function is specific to the Positive Poisson
#Intercept term moments: identical for all times
for(t in 1:T)
{
	gEst_i[count] = (mu_i[t]*(1-mu_i[t]*exp(-exp(eta_i[t]))))*y_ind_i[t]*(yvec_i[t]-mu_i[t])
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
		gEst_i[count] = (mu_i[t]*(1-mu_i[t]*exp(-exp(eta_i[t]))))*Zmat[subjectIndex+t-1,k]*y_ind_i[t]*(yvec_i[t]-mu_i[t])
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
				gEst_i[count] = (mu_i[s]*(1-mu_i[s]*exp(-exp(eta_i[s]))))*Xmat[subjectIndex+s-1,k]*y_ind_i[t]*(yvec_i[t]-mu_i[t])
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

} # end validMCPP_EC #


