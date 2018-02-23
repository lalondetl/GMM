

#' Multiple Testing p-Value Adjustment
#'
#' This function adjusts a list of p-values for multiple testing.  
#' @param pvals The vector of p-values to be adjusted.  
#' @keywords Multiple Testing
#' @export
#' @examples
#' pAdjusted()

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



