

#' Generalized Estimating Equations, Truncated Count component of Longitudinal Hurdle Model
#'
#' This function calculates the Generalized Estimating Equations (GEE) parameter estimates and standard errors for the truncated count component of a hurdle model for longitudinal excess zero count responses with independent working correlation structure, based on Dobbie and Welsh (2001).  Within the estimation method, responses are filtered to include only the positive counts from the original data.  Data must be organized by subject, and an intercept term is assumed for both the logistic and truncated count components of the model.  The function outputs a list with parameter estimates betaHat and parameter covariance estimate covEst, along with estimated rates muHat and residuals r_c.  
#' @param y	The vector of response counts, ordered by subject, time within subject.  
#' @param subjectID	The vector of subject ID values for each response.  
#' @param  N The number of subjects.  
#' @param X_c The design matrix for all covariates in the truncated count component of the model, including an intercept.  
#' @keywords GEE
#' @export
#' @examples
#' HurdleGEE_c()



HurdleGEE_c = function(y,subjectID,N,X_c){

# INVERSE LINK FUNCTION #
g_inv = function(x){exp(x)}

# NORM #
norm = function(x){sqrt(t(x)%*%x)}


# MINIMIZE EE USING ITERATIVE METHOD OF LIANG / ZEGER / QAQISH #
betaHat = rep(0,ncol(X_c))
deltaBeta = rep(10,ncol(X_c))
epsilon = 0.0001

while(norm(deltaBeta) > epsilon)
{
	# INITIALIZE INDEX, VALUE #
	index = 1
	sumA = matrix(0,ncol(X_c),ncol(X_c))
	sumB = rep(0,ncol(X_c))

	# CONSTRUCT DELTABETA COMPONENTS BY SUBJECT #
	for(i in 1:N)
	{
		# UPDATE RESPONSE, PREDICTORS, INDEX #
		y_i = as.vector(y[which(subjectID == subjectID[index])])
		X_c_i = as.matrix(X_c[which(subjectID == subjectID[index]),])
		index = max(which(subjectID == subjectID[index]))+1

		# SYSTEMATIC COMPONENT #
		eta_i = as.vector(X_c_i %*% betaHat)

		# ESTIMATED VALUES #
		lambda_i = as.vector(g_inv(eta_i))
		mu_i = as.vector(lambda_i / (1 - exp(-lambda_i)))
		y_ind_i = ifelse(y_i>0,1,0)

		# RESIDUAL VECTOR #
		c_i = y_ind_i*(y_i - mu_i)

		# WORKING COVARIANCE STRUCTURE #
		V_i = diag(mu_i*(1-lambda_i+mu_i))

		# DERIVATIVE MATRIX #
		D_i = lambda_i*X_c_i

		# UPDATE VALUES #
		sumA = sumA + t(D_i) %*% solve(V_i) %*% D_i
		sumB = sumB + t(D_i) %*% solve(V_i) %*% c_i
	}

	# UPDATE BETAHAT #
	deltaBeta = solve(sumA) %*% sumB
	betaHat = betaHat + deltaBeta
}


# OBTAIN STANDARD ERRORS #
##	USE BETAHAT	##

index=1
sumJ = matrix(0,ncol(X_c),ncol(X_c))
sumK = matrix(0,ncol(X_c),ncol(X_c))

for(i in 1:N)
{
	# UPDATE RESPONSE, PREDICTORS, INDEX #
	y_i = as.vector(y[which(subjectID == subjectID[index])])
	X_c_i = as.matrix(X_c[which(subjectID == subjectID[index]),])
	index = max(which(subjectID == subjectID[index]))+1

	# SYSTEMATIC COMPONENT #
	eta_i = as.vector(X_c_i %*% betaHat)

	# ESTIMATED VALUES #
	lambda_i = as.vector(g_inv(eta_i))
	mu_i = as.vector(lambda_i / (1 - exp(-lambda_i)))
	y_ind_i = ifelse(y_i>0,1,0)

	# RESIDUAL VECTOR #
	c_i = y_ind_i*(y_i - mu_i)

	# WORKING COVARIANCE STRUCTURE #
	V_i = diag(mu_i*(1-lambda_i+mu_i))

	# DERIVATIVE MATRIX #
	D_i = lambda_i*X_c_i

	# UPDATE VALUES #
	sumJ = sumJ + t(D_i) %*% solve(V_i) %*% D_i
	sumK = sumK + t(D_i) %*% solve(V_i) %*% c_i %*% t(c_i) %*% t(solve(V_i)) %*% D_i
}

covEst = solve(sumJ) %*% sumK %*% solve(sumJ)



# OBTAIN RESIDUALS #
lambdaHat = as.vector(g_inv(X_c %*% betaHat))
muHat = as.vector(lambdaHat / (1 - exp(-lambdaHat)))
r_c = ifelse(y>0,as.vector((y - muHat) / (sqrt(muHat*(1-lambdaHat+muHat)))),NA)

list(betaHat=betaHat, muHat=muHat, r_c=r_c,covEst=covEst)

} # END HurdleGEE_c #







