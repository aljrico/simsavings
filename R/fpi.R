#' fpi
#'
#'
#'
#' This function computes the time dependent value of pi, based on the present situation of the strategy.
#'
#' @export
#' @examples
#' print(fpi(A = 0.5, K = 42, X = 100, C = 20, titme = 33))
fpi <- function(A, K, X, C, time){
	g <- sum(C[-c(1:time)])
	xpi <- A*(K + X + g)
	return(xpi)
}
