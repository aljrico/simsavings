#' Compute Return
#'
#'
#'
#' This function computes the approximated yearly return from a time series of the accumulated wealth of a saving scheme.
#'
#' @param c The annual quantity that is added to the investment fund.
#' @param years The years the the savings scheme lasts. Inclulding both saving up and spending times.
#' @export
#' @examples
#' v <- c()
#' v[1] <- 10
#' for(i in 1:59) v[i+1] <- v[i]*1.01 + 10
#' compute_return(vec =  v, a = 0.05, years = 60)
compute_return <- function(vec, c = 10, years = 60){
	x_m <- median(vec)
	ret <- (1/years)*(-1 + (1 + (8*(x_m))/(c*years))^(1/2))*100
	return(ret)
}
