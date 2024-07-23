#' Compute RMSE with estimated beta
#'
#' @param beta  length(lambda) x p estimated beta matrix
#' @param  truebeta  true beta vector
#' @param  lambda  grid of lambda's
#'
#' @return RMSE value for each lambda
#' @export
#'
#' @examples
#' upper = 2; lower = -6
#' lambda =  2^seq(upper, lower, by = -0.2)
#' p = 50; len = length(lambda)
#' estimated_beta = matrix(rnorm(p*len, 0, 1), ncol = p)
#' true_beta = rnorm(50, 0, 1)
#' lasso.warm.rmse = rmse(estimated_beta, true_beta, lambda)
rmse = function(beta, truebeta, lambda) {
  rmse = matrix(NA, nrow=length(lambda), ncol=1)
  for(i in 1:length(lambda)){
    rmse[i,] = sqrt(sum((beta[i,]-truebeta)^2)/sum(truebeta^2))
  }
  return(rmse)
}
