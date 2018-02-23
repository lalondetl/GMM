

#' Two-Step Generalized Method of Moments, Longitudinal Bernoulli Outcome
#'
#' This function calculates the Generalized Method of Moments (GMM) parameter estimates and standard errors for longitudinal Bernoulli (0/1) responses.  The function allows for unbalanced data, meaning subjects can have different numbers of times of observation.  Both time-independent covariates and time-dependent covariates can be accommodated.  Time-dependent covariates can be handled either by specifying the type of each time-dependent covariate, or by allowing the data to determine appropriate moment conditions through the extended classification method.  
#' @param yvec The vector of responses, ordered by subject, time within subject.
#' @param subjectID The vector of subject ID values for each response.  
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param Tvec The vector of times for each subject.  
#' @param N The number of subjects.  
#' @param mc The method of identifying appropriate moment conditions, either 'EC' for extended classification (default) or 'Types' for user-identified types.  
#' @param covTypeVec The vector indicating the type of each time-dependent covariate, according to the order of the columns of Xmat.  
#' @keywords GMM
#' @export
#' @examples
#' TSGMM_Ber()




TSGMM_Ber = function(yvec,subjectID,Zmat,Xmat,Tvec,N,mc='EC',covTypeVec=c(-1)){

library("MASS")

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



######################################
# OBTAIN INITIAL PARAMETER ESTIMATES #
######################################

# INDEPENDENT GEE #
library("gee")

if(K0==0){ZX = Xmat}
if(K0!=0){ZX = cbind(Zmat,Xmat)}

betaI = gee(yvec ~ ZX,id=subjectID,family=binomial,corstr="independence")$coefficients

# ALL ZEROS #
# betaI = c(rep(0,K))

######################################
######################################


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
	types = validComb_EC_Ber(yvec,Zmat,Xmat,betaI,Tvec,alpha)
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

	# FILL G AND VN USING VALIDMCBER FOR EACH SUBJECT #
	for(i in 1:N)
	{
		subjectIndex = sum(Tvec[0:(i-1)])+1
		if(mc=='EC'){Est_i = validMCBer_EC(yvec,subjectIndex,Zmat,Xmat,beta,Tvec[i],Tmax,Count,types)}
		if(mc=='Types'){Est_i = validMCBer_Types(yvec,subjectIndex,Zmat,Xmat,covTypeVec,beta,Tvec[i],Tmax,Count)}

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
	W = ginv(VN / D)

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

	if(mc=='EC'){Est_i = validMCBer_EC(yvec,subjectIndex,Zmat,Xmat,betahat,Tvec[i],Tmax,Count,types)}
	if(mc=='Types'){Est_i = validMCBer_Types(yvec,subjectIndex,Zmat,Xmat,covTypeVec,betahat,Tvec[i],Tmax,Count)}

	gEst_i = Est_i[[1]]
	Count = Est_i[[2]]
	VN = VN + gEst_i%*%t(gEst_i)

	if(mc=='EC'){dBetagEst_i = validMDBer_EC(yvec,subjectIndex,Zmat,Xmat,betahat,Tvec[i],Tmax,types)}
	if(mc=='Types'){dBetagEst_i = validMDBer_Types(yvec,subjectIndex,Zmat,Xmat,covTypeVec,betahat,Tvec[i],Tmax)}
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
W = ginv(VN / D)

AsymptoticWeight = t(dBetaG) %*% W %*% dBetaG
AsymptoticCovariance = (1/N)*ginv(AsymptoticWeight)

########################################################################
########################################################################


list(betaHat=betahat, covEst = AsymptoticCovariance)

} # END TSGMM_BER #










###############################################################################
#validComb_EC_Ber function definition
#	Calculates the valid moment combinations using Extended 
#		Classification.  Allows unbalanced data.
#
#Input:	yvec		vector of binary responses.  organized by subject, by time within subject
#	Zmat		design matrix for time-independent covariates ((N*T) x K0)
#	Xmat		design matrix for time-dependent covariates ((N*T) x Ktv)
#	betaI		estimates of the model parameters (1+K0+Ktv x 1)
#	Tvec		vector of number of time points for each subject (Nx1)
#	alpha		significance level for inclusion
#
#Output:	types	Tmax x Tmax(Ktv) matrix, 1's indicating valid MC
###############################################################################
validComb_EC_Ber = function(yvec,Zmat,Xmat,betaI,Tvec,alpha){

####################
# DEFINE CONSTANTS #
####################

if(!is.matrix(Zmat)){K0 = 0}
else if(is.matrix(Zmat)){K0 = ncol(Zmat)}
Ktv = ncol(Xmat)
K = 1+K0+Ktv # TOTAL NUMBER OF PARAMETERS #

Tmax = max(Tvec)
N = length(Tvec)

TimePoint = rep(0,length(yvec))
subjectIndex = 1
for(i in 1:N)
{
	subjectIndex = sum(Tvec[0:(i-1)])+1
	TimePoint[(subjectIndex:(subjectIndex+Tvec[i]-1))] = seq(1:Tvec[i])
}

####################
####################


############################
# CALCULATE INITIAL VALUES #
############################

# MEAN AND SYSTEMATIC ESTIMATES #
# NOTE: THIS IS SPECIFIC TO THE BERNOULLI #

if(K0!=0){ZX = cbind(rep(1,nrow(Xmat)),Zmat,Xmat)}
if(K0==0){ZX = cbind(rep(1,nrow(Xmat)),Xmat)}

eta = ZX %*% betaI
mu = exp(eta)/(1+exp(eta))

# RESIDUALS #
r_raw = yvec - mu
r = vector(mode="list",length=Tmax)
for(t in 1:Tmax)
{
	r[[t]] = (r_raw[which(TimePoint==t)]-mean(r_raw[which(TimePoint==t)],na.rm=TRUE))/sqrt(var(r_raw[which(TimePoint==t)],na.rm=TRUE))
}

##################################
##################################

types = matrix(0,Tmax,Tmax*Ktv)

################################################
# CALCULATE CORRELATION FOR EACH TDC PARAMETER #
################################################

for(j in 1:Ktv)
{
	# FIND DERIVATIVES #
	dBetamu_j_raw = mu*(1-mu)*Xmat[,j] 
	d_j = vector(mode="list",length=Tmax)

	for(t in 1:Tmax)
	{
		d_j[[t]] = (dBetamu_j_raw[which(TimePoint==t)]-mean(dBetamu_j_raw[which(TimePoint==t)],na.rm=TRUE))/sqrt(var(dBetamu_j_raw[which(TimePoint==t)],na.rm=TRUE))
	}

	# CALCULATE CORRELATION AND MIXED MOMENT #
	rho_j = matrix(0,Tmax,Tmax)
	mu22_j = matrix(0,Tmax,Tmax)
	for(t in 1:Tmax){for(s in 1:Tmax){
		if(s<t)
		{	
			rho_j[s,t] = cor(d_j[[s]][which(Tvec[which(Tvec>=s)]>=t)],r[[t]],use = "na.or.complete")
			mu22_j[s,t] = (1/N)*sum((d_j[[s]][which(Tvec[which(Tvec>=s)]>=t)])^2*(r[[t]])^2,na.rm=TRUE)
		}
		
		if(s>=t)
		{	
			rho_j[s,t] = cor(d_j[[s]],r[[t]][which(Tvec[which(Tvec>=t)]>=s)],use = "na.or.complete")
			mu22_j[s,t] = (1/N)*sum((d_j[[s]])^2*(r[[t]][which(Tvec[which(Tvec>=t)]>=s)])^2,na.rm=TRUE)
		}
	}}

	rhoStar_j = rho_j/(sqrt(mu22_j/N))

	# SELECT NON-SIGNIFICANT PAIRS FOR EACH PARAMETER #
	types_j = matrix(0,Tmax,Tmax)
	# types_j[abs(rhoStar_j) < qnorm(1-(alpha/2))] = 1

	# FIND P-VALUES, APPLY P_ACT #
	pvals = as.vector(2*pnorm(abs(rhoStar_j),lower.tail=FALSE))
	adj_pvals = matrix(c(pAdjusted(pvals)),Tmax,Tmax)

	types_j[adj_pvals>=0.05]=1

	types[,((Tmax*(j-1)+1):(Tmax*(j-1)+Tmax))] = types_j
}

################################################
################################################

types

} # END validComb_EC_Ber #






###############################################################################
#pACT function definition
#	Adjusts the first of a list of p-values for multiple testing 
#
#Input:		pvals		a list of p-values
#
#Output:	p_ACT		an adjusted *minimum* p-value
###############################################################################
pACT = function(pvals){

library(mvtnorm)
level=25000

L = length(pvals)

# FOR THE MINIMUM.  REPEATED WITHOUT PMIN ITERATIVELY #

minp=min(pvals)

if (minp==0) {p_ACT=0}
if (minp>=.5) {p_ACT=1}

if (minp>0 & minp < .5)
{
	lower=rep(qnorm(minp/2),L)
	upper=rep(qnorm(1-minp/2),L)

	# CALCULATE COVARIANCE MATRIX v #
	v = diag(1,length(pvals))

	p_ACT=1-pmvnorm(lower=lower,upper=upper,sigma=v,maxpts=level,abseps=.0000000000001)
}

p_ACT
} # END pACT #





###############################################################################
#pAdjusted function definition
#	Adjusts a list of p-values for multiple testing 
#
#Input:		pvals		a list of p-values
#
#Output:	p_Adjusted	an adjusted list of p-values
###############################################################################
pAdjusted = function(pvals){

p_Adjusted = rep(0,length(pvals))
order_index = order(pvals)

new_pvals = pvals

for(i in 1:length(pvals))
{
	p_Adjusted[order_index[i]] = pACT(new_pvals)
	new_pvals = new_pvals[-order(new_pvals)[1]]
}

p_Adjusted

}









###############################################################################
#validMCBer_EC function definition
#	Calculates the valid moment conditions for one subject.  Allows unbalanced data.
#
#Input:	
#	yvec		vector of responses.  organized by subject, by time within subject
#	subjectIndex	location subject i within yvec
#	Zmat		design matrix for time-independent covariates ((N*T) x K0)
#	Xmat		design matrix for time-dependent covariates ((N*T) x Ktv)
#	betaI		estimates of the model parameters (1+K0+Ktv x 1)
#	T		number of time points, subject i
#	Tmax		maximum number of time points observed
#	Count		vector of running count of valid moment conditions for all subjects
#	types		indicates valid combinations for each TDC parameter

#Output:
#	g_est_i		vector of valid moment conditions for i^th subject
#	Count		incremented for all valid moment conditions of subject i
###############################################################################
validMCBer_EC = function(yvec,subjectIndex,Zmat,Xmat,betaI,T,Tmax,Count,types){

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
	gEst_i[count] = (mu_i[t]/(1+exp(eta_i[t])))*(yvec_i[t]-mu_i[t])
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
		gEst_i[count] = (mu_i[t]/(1+exp(eta_i[t])))*Zmat[subjectIndex+t-1,k]*(yvec_i[t]-mu_i[t])
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
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(yvec_i[t]-mu_i[t])
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

} # end validMCBer_EC #










###############################################################################
#validMCBer_Types function definition
#	Calculates the valid moment conditions for one subject.  Allows unbalanced data.
#
#Input:	
#	yvec		vector of responses.  organized by subject, by time within subject
#	subjectIndex	location subject i within yvec
#	Zmat		design matrix for time-independent covariates ((N*T) x K0)
#	Xmat		design matrix for time-dependent covariates ((N*T) x Ktv)
#	covTypeVec	vector indicating the type of each time-dependent covariate
#	betaI		estimates of the model parameters (1+K0+Ktv x 1)
#	T		number of time points, subject i
#	Tmax		maximum number of time points observed
#	Count		vector of running count of valid moment conditions for all subjects

#Output:
#	g_est_i		vector of valid moment conditions for i^th subject
#	Count		incremented for all valid moment conditions of subject i
###############################################################################
validMCBer_Types = function(yvec,subjectIndex,Zmat,Xmat,covTypeVec,betaI,T,Tmax,Count){

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




##########################################################
# DEFINE VECTOR OF VALID MOMENT CONDITIONS FOR SUBJECT i #
##########################################################

gEst_i = rep(0,Lmax)
count = 1

#NOTE: gEst_i function is specific to the Bernoulli / Binomial
#Intercept term moments: identical for all times
for(t in 1:T)
{
	gEst_i[count] = (mu_i[t]/(1+exp(eta_i[t])))*(yvec_i[t]-mu_i[t])
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
		gEst_i[count] = (mu_i[t]/(1+exp(eta_i[t])))*Zmat[subjectIndex+t-1,k]*(yvec_i[t]-mu_i[t])
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
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(yvec_i[t]-mu_i[t])
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
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(yvec_i[t]-mu_i[t])
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
			gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(yvec_i[s]-mu_i[s])
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
				gEst_i[count] = (mu_i[s]/(1+exp(eta_i[s])))*Xmat[subjectIndex+s-1,k]*(yvec_i[t]-mu_i[t])
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

} # end validMCBer_Types #











###############################################################################
#validMDBer_EC function definition
#	Calculates the derivative of the valid moment conditions vector for one subject.  Allows unbalanced data.
#
#Input:	
#	yvec		vector of responses.  organized by subject, by time within subject
#	subjectIndex	location of subject i within yvec
#	Zmat		design matrix for time-independent covariates ((N*T) x K0)
#	Xmat		design matrix for time-dependent covariates ((N*T) x Ktv)
#	betaI		estimates of the model parameters (1+K0+Ktv x 1)
#	N		number of subjects
#	T		number of time points, subject i
#	Tmax		maximum number of time points observed
#	types		matrix indicating valid MC's based on Extended Classification

#Output:
#	dBetagEst_i	derivative of vector of valid moment conditions for i^th subject
###############################################################################
validMDBer_EC = function(yvec,subjectIndex,Zmat,Xmat,betaI,T,Tmax,types){

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
	types_j = types[,((Tmax*(j-1)+1):(Tmax*(j-1)+Tmax))]

	for(s in 1:T){
		for(t in 1:T){
			if(types_j[s,t] == 1){
				for(k in 1:K)
				{
					dBetag_i[count,k] = (-1)*dBetamu_i[s,1+K0+j]*dBetamu_i[t,k] + (Xmat[(subjectIndex+s-1),j])*d2Betamu_i_part[s,k]*(yvec_i[t]-mu_i[t])
				}
				count = count + 1
			}
		}
	
		# UPDATE COUNT FOR MISSING TIMES #
		if(T<Tmax){count = count + sum(types_j[s,((T+1):Tmax)])}
	}

	# UPDATE COUNT FOR MISSING TIMES #
	if(T<Tmax){count = count + sum(types_j[((T+1):Tmax),])}

} # END TIME-DEPENDENT COVARIATE DERIVATIVES


# dBetag_i IS NOW A (L x K) MATRIX OF DERIVATIVES FOR SUBJECT i #

dBetag_i
} # END validMDBer_EC #








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






