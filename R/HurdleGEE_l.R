

#' Generalized Estimating Equations, Logistic component of Longitudinal Hurdle Model
#'
#' This function calculates the Generalized Estimating Equations (GEE) parameter estimates and standard errors for the logistic component of a hurdle model for longitudinal excess zero count responses with independent working correlation structure, based on Dobbie and Welsh (2001).  Responses are treated as binary indicators of 0 count values, therefore the probability modeled is that of a zero count.  Data must be organized by subject, and an intercept term is assumed for both the logistic and truncated count components of the model.  The function outputs a list with parameter estimates betaHat and parameter covariance estimate covEst, along with estimated probabilities pHat and residuals r_l.  
#' @param y	The vector of response counts, ordered by subject, time within subject.  
#' @param subjectID	The vector of subject ID values for each response.  
#' @param  N The number of subjects.  
#' @param X_l The design matrix for all covariates in the logistic component of the model, including an intercept.  
#' @keywords GEE
#' @export
#' @examples
#' HurdleGEE_l()



HurdleGEE_l = function(y,subjectID,N,X_l){

library(MASS) 

# CREATE BINARY RESPONSE #
# MODEL CERTAIN ZERO!!! #
y_binary = ifelse(y>0,0,1)

# INVERSE LINK FUNCTION #
g_inv = function(x){exp(x)/(1+exp(x))}

# NORM #
norm = function(x){sqrt(t(x)%*%x)}


# MINIMIZE EE USING ITERATIVE METHOD OF LIANG / ZEGER / QAQISH #
betaHat = rep(0,ncol(X_l))
deltaBeta = rep(10,ncol(X_l))
epsilon = 0.0001

while(norm(deltaBeta) > epsilon)
{
	# INITIALIZE INDEX, VALUE #
	index = 1
	sumA = matrix(0,ncol(X_l),ncol(X_l))
	sumB = rep(0,ncol(X_l))

	# CONSTRUCT DELTABETA COMPONENTS BY SUBJECT #
	for(i in 1:N)
	{
		# UPDATE RESPONSE, PREDICTORS, INDEX #
		y_binary_i = as.vector(y_binary[which(subjectID == subjectID[index])])
		X_l_i = as.matrix(X_l[which(subjectID == subjectID[index]),])
		index = max(which(subjectID == subjectID[index]))+1

		# SYSTEMATIC COMPONENT #
		eta_i = as.vector(X_l_i %*% betaHat)

		# ESTIMATED VALUES #
		p_i = as.vector(g_inv(eta_i))

		# RESIDUAL VECTOR #
		b_i = y_binary_i - p_i

		# WORKING COVARIANCE STRUCTURE #
		V_i = diag(p_i*(1-p_i))

		# DERIVATIVE MATRIX #
		D_i = p_i*(1-p_i)*X_l_i

		# UPDATE VALUES #
		sumA = sumA + t(D_i) %*% solve(V_i) %*% D_i
		sumB = sumB + t(D_i) %*% solve(V_i) %*% b_i
	}

	# UPDATE BETAHAT #
	deltaBeta = solve(sumA) %*% sumB
	betaHat = betaHat + deltaBeta
}


# OBTAIN STANDARD ERRORS #
##	USE BETAHAT	##

index=1
sumJ = matrix(0,ncol(X_l),ncol(X_l))
sumK = matrix(0,ncol(X_l),ncol(X_l))

for(i in 1:N)
{
	# UPDATE RESPONSE, PREDICTORS, INDEX #
	y_binary_i = as.vector(y_binary[which(subjectID == subjectID[index])])
	X_l_i = as.matrix(X_l[which(subjectID == subjectID[index]),])
	index = max(which(subjectID == subjectID[index]))+1

	# SYSTEMATIC COMPONENT #
	eta_i = as.vector(X_l_i %*% betaHat)

	# ESTIMATED VALUES #
	p_i = as.vector(g_inv(eta_i))

	# RESIDUAL VECTOR #
	b_i = y_binary_i - p_i

	# WORKING COVARIANCE STRUCTURE #
	V_i = diag(p_i*(1-p_i))

	# DERIVATIVE MATRIX #
	D_i = p_i*(1-p_i)*X_l_i

	# UPDATE VALUES #
	sumJ = sumJ + t(D_i) %*% solve(V_i) %*% D_i
	sumK = sumK + t(D_i) %*% solve(V_i) %*% b_i %*% t(b_i) %*% t(solve(V_i)) %*% D_i
}

covEst = solve(sumJ) %*% sumK %*% solve(sumJ)


# OBTAIN RESIDUALS #
pHat = as.vector(g_inv(X_l %*% betaHat))
r_l = as.vector((y_binary - pHat) / (sqrt(pHat*(1-pHat))))

list(betaHat=betaHat, pHat=pHat, r_l=r_l,covEst=covEst)

} # END HurdleGEE_l #






