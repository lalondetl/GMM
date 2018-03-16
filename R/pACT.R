

#' Actual Adjusted Minimum p-Value
#'
#' This function takes a list of ordered p-values and returns an adjusted minimum p-value.  
#' @param pvals The vector of p-values to be adjusted.  
#' @keywords Multiple Testing
#' @export
#' @examples
#' pACT()

pACT = function(pvals){

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

	p_ACT=1-mvtnorm::pmvnorm(lower=lower,upper=upper,sigma=v,maxpts=level,abseps=.0000000000001)
}

p_ACT
} # END pACT #





