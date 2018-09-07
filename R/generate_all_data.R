#' Generate all data
#'
#'
#'
#' This function generates a data frame containing the results from all the simulations this package is capable of..
#'
#' @param alpha Expected return of the risky market
#' @param sigma Standard deviation of the risky market
#' @param years Duration of the savings scheme.
#' @param include.mortality Boolean that decides whether the simulation contains mortality.
#' @export
#' @examples
#' head(generate_all_data())

generate_all_data <- function(
	alpha = 0.343,
	sigma = 0.1544,
	a = 10,
	years = 60,
	nsim = 1e3,
	pi = 0.1,
	K = 42,
	A = 0.5,
	include.mortality = FALSE
){

	# CPPI simple --------------------------------------------------------------------
	cat("... CPPI Simple ... \n")
	X_T <- cppi_c(pi = pi,
								nsim = nsim,
								alpha = alpha,
								sigma = sigma,
								a = a,
								years = years)

	final_wealth <- as_tibble(as.data.frame(X_T))
	final_wealth$model <- "cppi-simple"
	final_wealth$ret <- (1/60)*(-1 + (1 + (8*(X_T))/(a*60))^(1/2))*100

	write.csv(X_T, "data/cppi_simple.csv")

	K = ES(X_T)*factor_kes(A =A)


	# Alternative simple ------------------------------------------------------
	cat("... Alternative Simple ... \n")
	X_T <- alt_c(K = K,
							 nsim = nsim,
							 alpha = alpha,
							 sigma = sigma,
							 a = a,
							 years = years,
							 A_factor = A)

	df <- as_tibble(as.data.frame(X_T))
	df$model <- "alt-simple"
	df$ret <- (1/60)*(-1 + (1 + (8*(X_T))/(a*60))^(1/2))*100
	final_wealth <- rbind(final_wealth,df)

	write.csv(X_T, "data/alt_simple.csv")

	if(include.mortality == TRUE){
		# CPPI | Mortality --------------------------------------------------------
		cat("... CPPI Mortality ... \n")
		X_T <- cppi_mort_fasto(pi = pi,
													 nsim = nsim,
													 alpha = alpha,
													 sigma = sigma,
													 a = a,
													 years = years)
		df <- as_tibble(as.data.frame(X_T))
		df$model <- "cppi-mort"
		df$ret <- (1/60)*(-1 + (1 + (8*(X_T))/(a*60))^(1/2))*100
		final_wealth <- rbind(final_wealth,df)

		write.csv(X_T, "data/cppi_mortality.csv")

		# K = ES(X_T)*factor_kes(A = A)
		K = es_to_k(A = A, pi = pi, nsim = 1e2, k_max = 1000, err = 0.01)
		if(is.na(K)) K <- es_to_k(A = A, pi = pi, nsim = 1e2, k_max = 1000, size = 10, err = 0.1)
		if(is.na(K)) K <- es_to_k(A = A, pi = pi, nsim = 1e2, k_max = 10000, size = 25, err = 0.25)
		if(is.na(K)) K <- es_to_k(A = A, pi = pi, nsim = 1e2, k_max = 1000, size = 50, err = 0.5)

		# Alternative | Mortality ------------------------------------------------------
		cat("... Alternative Mortality ... \n")
		X_T <- alt_mort_fasto(K = K,
													nsim = nsim,
													alpha = alpha,
													sigma = sigma,
													a = a,
													years = years)
		df <- as_tibble(as.data.frame(X_T))
		df$model <- "alt-mort"
		df$ret <- (1/60)*(-1 + (1 + (8*(X_T))/(a*60))^(1/2))*100
		final_wealth <- rbind(final_wealth,df)

		write.csv(X_T, "data/alt_mortality.csv")

	}

	return(final_wealth)
}
