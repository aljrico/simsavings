#' Expected Shortfall (ES)
#'
#'
#'
#' This function computes the Expected Shortfall, provided a wealth distribution.
#'
#' @param a The chosen quantile of the distribution where the Value at Risk (VaR) is set.
#' @export
#' @examples
#' ES(distr = rnorm(1000), a = 0.05)
ES <- function(distr, a = 0.05){
	VaR <- quantile(distr, a, na.rm = TRUE)
	ES <- mean(distr[distr<VaR])

	return(ES)
}
