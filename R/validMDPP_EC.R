

#' Generalized Method of Moments Valid Moment Combination Derivatives for one Subject of Longitudinal Zero-Truncated Count Responses, using Extended Classification
#' 
#' This function calculates the values of the derivatives of all valid moment conditions for a single subject in a longitudinal study with zero-truncated count ("positive Poisson") outcomes.  It allows for unbalanced data, and uses the extended classification method to determine validity of moment conditions.  The function returns a matrix of derivatives for all valid moment condition for subject i.  
#' @param yvec The vector of responses, ordered by subject, time within subject.
#' @param subjectIndex The location of the first index of subject i responses within yvec.  
#' @param Zmat The design matrix for time-independent covariates ((N*T) x K0).  
#' @param Xmat The design matrix for time-dependent covariates ((N*T) x Ktv).  
#' @param betaI The current or initial estimates of the model parameters (1+K0+Ktv x 1).  
#' @param N The number of subjects.  
#' @param T The number of time points for subject i.  
#' @param Tmax The maximum number of times of observation among all subjects.  
#' @param types A matrix indicating valid moment conditions for all time-dependent covariate parameters.  
#' @keywords GMM
#' @export
#' @examples
#' validMDPP_EC()



validMDPP_EC = function(yvec,subjectIndex,Zmat,Xmat,betaI,T,Tmax,types){

####################
# DEFINE CONSTANTS #
####################

if(!is.matrix(Zmat)){K0 = 0}
else if(is.matrix(Zmat)){K0 = ncol(Zmat)}
Ktv = ncol(Xmat)
K = 1+K0+Ktv # TOTAL NUMBER OF PARAMETERS #

# MAXIMUM NUMBER OF ESTIMATING EQUATIONS (PER SUBJECT) #
Lmax = 1*Tmax + K0*Tmax  + sum(types)

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







#########################################################
#	DEFINE DERIVATIVE MATRICES FOR SUBJECT i	#
#		UNIQUE TO POSITIVE POISSON!!		##########################################################


##	FIRST DERIVATIVE OF MEAN	##

dBetamu_i = matrix(0,T,K)
for(t in 1:T)
{
dCount = 1

	# INTERCEPT #
	dBetamu_i[t,dCount] = (1)*(mu_i[t]*(1-mu_i[t]*exp(-exp(eta_i[t]))))
	dCount = dCount+1

	# TIC #
	if(K0!=0)
	{
	for(j in 1:K0)
	{
		dBetamu_i[t,dCount] = (Zmat[subjectIndex+t-1,j])*(mu_i[t]*(1-mu_i[t]*exp(-exp(eta_i[t]))))
		dCount = dCount+1
	}
	}

	# TDC #
	for(j in 1:ncol(Xmat))
	{
		dBetamu_i[t,dCount] = (Xmat[subjectIndex+t-1,j])*(mu_i[t]*(1-mu_i[t]*exp(-exp(eta_i[t]))))
		dCount = dCount+1
	}
}

##	PART OF SECOND DERIVATIVE OF MEAN	##
#		(Xmat[i,t,j] OMITTED)		##

d2Betamu_i_part = matrix(0,T,K)
for (t in 1:T)
{
	for (k in 1:K)
	{
		d2Betamu_i_part[t,k] = (2*mu_i[t]*exp(-exp(eta_i[t]))*exp(eta_i[t])*Xmat[subjectIndex+t-1,k]) + (dBetamu_i[t,k]*(1-2*mu_i[t]*exp(-exp(eta_i[t]))))
	}
}




#################################################
#	MATRIX OF DERIVATIVES FOR SUBJECT i	#
#		INDICATOR FOR PP!!		##################################################

dBetag_i = matrix(0,Lmax,K)
count = 1

## INTERCEPT EE'S: IDENTICAL FOR ALL TIMES ##
for(t in 1:T)
{
	s = t
	j = 1
	for(k in 1:K)
	{
		dBetag_i[count,k] = (-1)*dBetamu_i[s,j]*dBetamu_i[t,k]*y_ind_i[t] + d2Betamu_i_part[s,k]*(1)*y_ind_i[t]*(yvec_i[t]-mu_i[t])
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
			dBetag_i[count,k] = (-1)*dBetamu_i[s,1+j]*dBetamu_i[t,k]*y_ind_i[t] + (Zmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*y_ind_i[t]*(yvec_i[t]-mu_i[t])
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
	types_j = types[,((Tmax*(j-1)+1):(Tmax*(j-1)+Tmax))]

	for(s in 1:T){
		for(t in 1:T){
			if(types_j[s,t] == 1){
				for(k in 1:K)
				{
					dBetag_i[count,k] = (-1)*dBetamu_i[s,1+K0+j]*dBetamu_i[t,k]*y_ind_i[t] + (Xmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*y_ind_i[t]*(yvec_i[t]-mu_i[t])
				}
				count = count + 1
			}
		}
	
		# UPDATE COUNT FOR MISSING TIMES #
		if(T<Tmax){count = count + sum(types_k[s,((T+1):Tmax)])}
	}

	# UPDATE COUNT FOR MISSING TIMES #
	if(T<Tmax){count = count + sum(types_k[((T+1):Tmax),])}

} # END TIME-DEPENDENT COVARIATE DERIVATIVES


# dBetag_i IS NOW A (L x K) MATRIX OF DERIVATIVES FOR SUBJECT i #

dBetag_i
} # END validMDPP_EC #



