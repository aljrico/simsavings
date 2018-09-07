#' From ES to K
#'
#'
#'
#' This function estimates the necessary value of K in order to make the Alternative Strategy equal the Expected Shortfall of the CPPI strategy.
#'
#' @param A Parameter A that loosely defines the risk aversion of the investor.
#' @param pi Constant proportion of the portfolio invested in risky assets (CPPI).
#' @param nsim Number of simulations.
#' @param err Maximum relative error allowed for the estimation.
#' @param k_max Upper boundary to the numeric method.
#' @param size = Width of the steps of the numerical method.
#' @export
#' @examples
#' es_to_k(nsim = 1e2)
es_to_k <- function(A = 0.5, pi = 0.7, nsim = 1e2, err = 0.01, k_max = 1000, size = 1){
	k_result <- c()
	es_init <- c()
	es_final <- c()
	es_record <- c()
	k_record <- seq(1, k_max, by = size)
	pi_record <- pi

	pb <- progress_bar$new(total = length(k_record))

	for(i in 1:length(pi_record)){
		es_init[i] <- cppi_mort_fasto(pi = pi_record[i], nsim =nsim) %>% na.omit() %>% ES()

		for(j in 1:length(k_record)) {
			pb$tick()
			es_record[j] <- alt_mort_fasto(K = k_record[j], nsim = nsim, A = A) %>% na.omit() %>% ES()
		}
		k_result[i] <- mean(which(abs((es_record - es_init[i])/(es_init[i])) < err))
	}
	return(k_result)
}
