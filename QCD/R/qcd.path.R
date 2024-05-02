#' Solve penalized quantile regression by pathwise coordinate descent for a grid of lambda's with user chosen penalty
#'
#' @param x n x p design matrix X
#' @param y n x 1 vector of response variables Y
#' @param tau quantile
#' @param funname user chosen "LASSO", "SCAD", and "MCP" as the penalty
#' @param a threshold parameter that adjusts constant penalty part for SCAD and MCP (a > 2 for SCAD, a > 1 for MCP)
#' @param weights can input different weight for each p coefficients (default is no weights)
#' @param lambda vector of lambda gird
#' @param nlambda number of lambda grid. automatically calculated if there is no input.
#' @param nudge whether nudge will be added when changing lambda (default is nudge = FALSE). recommend set.seed() before the function if nudge is used.
#' @param nudgesd user chosen standard deviation of the random nudge (default is 0.01)
#' @param standardize whether X will be scaled (default is standardize = TRUE)
#' @param thresh threshold of checking whether the coefficients converged
#' @param maxit maximum iteration for convergence
#' @param verbose whether the iteration number will be printed. verbose = TRUE will print the iteration.
#'
#' @return \item{beta}{a n x 1 matrix of coefficients, stored in sparse matrix format}
#' @return \item{dim}{dimension of coefficient matrix}
#' @return \item{lambda}{lambda value used}
#' @return \item{nobs}{number of observations}
#' @return \item{df}{number of nonzero coefficients}
#' @export
#'
#' @examples
#' n = 30; p = 30
#' x = array(rnorm(n*p), c(n,p))
#' for (j in 1:p){x[,j] = x[,j]/(norm(as.matrix(x[,j]), type="f")/sqrt(n))}
#' e = rnorm(n)
#' b = c(1, -1, rep(0, p-2))
#' y = x %*% b + e
#'
#' upper = 2; lower = -6
#' lambda =  2^seq(upper, lower, by = -0.2)
#'
#' qr.lasso.warm = qcd.path(x = x, y = y, tau = 0.5, funname = "LASSO",
#' weights = NULL, lambda = lambda, nlambda = length(lambda),
#' nudge = FALSE, standardize = TRUE,
#' thresh = 1e-06, maxit = 10000, verbose = TRUE)
#'
#' set.seed(1)
#' qr.lasso.warm.nudge = qcd.path(x = x, y = y, tau = 0.5, funname = "LASSO",
#' weights = NULL, lambda = lambda, nlambda = length(lambda),
#' nudge = TRUE, nudgesd = 0.01, standardize = TRUE,
#' thresh = 1e-06, maxit = 10000, verbose = TRUE)
qcd.path <- function(x, y, tau, funname = c("LASSO", "SCAD", "MCP"), a,
                     weights = NULL, lambda = NULL, nlambda = 30,
                     nudge = FALSE, nudgesd = 0.01, standardize = TRUE,
                     thresh = 1e-6, maxit = 100000, verbose = TRUE) {

  # store dimension
  np = dim(x)
  if (is.null(np) || (np[2] <= 1)) stop("x should be a matrix with 2 or more columns")
  n = as.integer(np[1]); p = as.integer(np[2])

  # get feature variable names
  vnames=colnames(x)
  if (is.null(vnames)) vnames = paste("V",seq(p),sep="")

  # check weights
  if (is.null(weights)) weights = rep(1,n)
  else if (length(weights) != n)
    stop(paste("Number of elements in weights (",length(weights),
               ") not equal to the number of rows of x (",n,")",sep=""))
  weights <- as.double(weights)

  # check standardize
  if (standardize) xs = scale(x)
  else xs = x

  # work out lambda values
  nlam = as.integer(nlambda)
  user_lambda = FALSE   # did user provide their own lambda values?
  if (is.null(lambda)) stop("lambda needs to be specified")
  else {  # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam = as.double(sort(lambda, decreasing = TRUE))
    nlam = as.integer(length(lambda))
  }

  # check if warm start will be used and create the path
  path = matrix(NA, nrow=length(lambda), ncol=p) # create matrix to store the path for different lambdas

  fit = NULL
  #set.seed(9)
  for (k in 1:nlam) {
    cur_lambda = ulam[k]
    cat("Fitting lambda index", k, ":", round(log2(cur_lambda), 2), fill = TRUE) # print which lambda we are at

    betahat = matrix(rep(NA, p*maxit), ncol=p)
    fit = qcd.fit(x = xs, y, tau, weights, cur_lambda, a,
                  funname = funname, warm = fit,
                  thresh = thresh, maxit = maxit, verbose = FALSE)
    path[k,] = fit$beta

    # If nudge = TRUE, we can choose the sd of rnorm() by nudgesd, recommend to set.seed before qcd.path()
    if (nudge == TRUE) fit = path[k,] + rnorm(p, mean = 0, sd = nudgesd)
    else fit = path[k,]
  }

  # output
  out = list(beta = path, dim = dim(path), lambda = ulam, nobs = n, df = rowSums(abs(path) >0))
  class(out) = "qcd.path"

  return(out)
}
