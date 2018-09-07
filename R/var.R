#' Value at Risk (VaR)
#'
#'
#'
#' This function computes the Value at Risk, provided a wealth distribution.
#'
#' @param a The chosen quantile of the distribution where the Value at Risk (VaR) is set.
#' @export
#' @examples
#' VaR(distr = rnorm(1000), a = 0.05)
VaR <- function(distr, a = 0.05){quantile(distr, a, na.rm = TRUE)}

