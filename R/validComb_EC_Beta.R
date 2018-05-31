

#' Generalized Method of Moments Valid Moment Combinations for Longitudinal Proportion Responses, using Extended Classification
#' 
#' This function calculates the values of valid moment combinations for two-step Generalized Method of Moments using the extended classification method, applied to longitudinal data with proportion outcomes.  It allows for unbalanced longitudinal data, meaning subjects can be observed for different numbers of times.  The function returns a vector "types" indicating validity of different moments conditions.  
#' @param yvec The vector of responses, ordered by subject, time within subject.
#' @param Zmat The design matrix for time-independent covariates.  
#' @param Xmat The design matrix for time-dependent covariates.  
#' @param betaI The current or initial estimates of the model parameters.  
#' @param Tvec The vector of times for each subject. 
#' @param alpha The significance level for inclusion of moment conditions.  
#' @keywords GMM
#' @export
#' @examples
#' validComb_EC_Beta()

validComb_EC_Beta = function(yvec,Zmat,Xmat,betaI,Tvec,alpha){

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
# NOTE: THIS IS SPECIFIC TO THE BETA #

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

} # END validComb_EC_Beta #



