#' KES Factor (factor_kes)
#'
#'
#'
#' This function computes the factor that theoretically relates K to Expected Shortfall.
#'
#' @param alpha Expected return of the risky market
#' @param sigma Standard deviation of the risky market
#' @param years Duration of the savings scheme.
#' @param theta 1 - the quantile of the Expected Shortfall
#' @export
#' @examples
#' factor_kes(A = 0.5)
factor_kes <- function(alpha = 0.0343,
											 sigma = 0.1544,
											 years = 60,
											 theta = 0.95,
											 A){
	factor <- 1/(-1 + (1/(1 - theta))*exp(alpha*A*years)*pnorm(qnorm(1-theta)- A*sigma*sqrt(years)))
}
