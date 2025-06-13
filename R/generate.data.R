#' Generate Peng and Wang (2015) simulation data
#'
#' @param n number of observations.
#' @param p number of variables. (p > 20)
#' @param signal signal of true beta's.
#' @param tau quantile value between 0 and 1.
#'
#' @return \item{Y}{n x 1 vector.}
#' @return \item{X}{n x p matrix.}
#' @return \item{true_beta}{p x 1 true beta vector.}
#' @import mvtnorm
#' @export
#'
#' @examples
#' set.seed(1)
#' data = generate.data(n = 150, p = 150, signal = 1, tau = 0.5)
#'
generate.data <- function(n, p, signal = 1, tau = 0.5) {
  indx <- c(6, 12, 15, 20)
  true_beta <- rep(0, p)
  true_beta[indx] <- signal
  true_beta[1] <- 0.7*qnorm(tau, mean = 0, sd = 1)
  Sigma <- 0.5^abs(outer(1:p, 1:p, '-'))
  X <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
  epsilon <- rnorm(n)
  Y <- X[,6] + X[,12] + X[,15] + X[,20] + 0.7 * pnorm(X[,1]) * epsilon
  return(list(Y = Y, X = X, true_beta = true_beta))
}


