#' Solve penalized quantile regression by coordinate descent for a single value of lambda with user chosen penalty
#'
#' @param x n x p design matrix X
#' @param y n x 1 vector of response variables Y
#' @param tau quantile
#' @param weights can input different weight for each p coefficients (default is no weights)
#' @param lambda single lambda value
#' @param a threshold parameter that adjusts constant penalty part for SCAD and MCP (a > 2 for SCAD, a > 1 for MCP)
#' @param funname user chosen "LASSO", "SCAD", and "MCP" as the penalty
#' @param warm whether warm-start will be used (default is NULL). do not use this argument when only running qcd.fit(). it is used when constructing regularization path with qcd.path() later.
#' @param thresh threshold of checking whether the coefficients converged
#' @param  maxit  maximum iteration for convergence (default is 100000)
#' @param  verbose  whether the iteration number will be printed. verbose = TRUE will print the iteration.
#'
#' @return \item{beta}{a n x 1 matrix of coefficients, stored in sparse matrix format}
#' @return \item{dim}{dimension of coefficient matrix}
#' @return \item{lambda}{lambda value used}
#' @return \item{df}{number of nonzero coefficients}
#' @import utils
#' @import stats
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
#' qr.lasso.choice = qcd.fit(x = x, y = y, tau = 0.5, lambda = 0.8,
#' funname = "LASSO", weights = NULL, warm = NULL,
#' thresh = 1e-06, maxit = 100, verbose = TRUE)

#' qr.scad.choice = qcd.fit(x = x, y = y, tau = 0.5, lambda = 0.8, a = 2.2,
#'                          funname = "SCAD", weights = NULL, warm = NULL,
#'                          thresh = 1e-06, maxit = 100, verbose = TRUE)

#' qr.mcp.choice = qcd.fit(x = x, y = y, tau = 0.5, lambda = 0.8, a = 2.2,
#'                         funname = "MCP", weights = NULL, warm = NULL,
#'                         thresh = 1e-06, maxit = 100, verbose = TRUE)
qcd.fit = function(x, y, tau, weights, lambda, a, funname = c("LASSO", "SCAD", "MCP"),
                   warm = NULL, thresh = 1e-6, maxit = 100000, verbose = FALSE) {

  if (is.null(warm)) fit = NULL
  else fit = warm # warm start input

  # run coordinate descent for different penalties
  if (funname == "LASSO") {
    fit = lasso.fit(x, y, tau, lambda, weights = weights, warm = fit,
                    thresh = thresh, maxit = maxit, verbose = verbose)
  }

  if (funname == "SCAD") { # total number of obs for SCAD : n+5
    fit = scad.fit(x, y, tau, lambda, a, weights = weights, warm = fit,
                   thresh = thresh, maxit = maxit, verbose = verbose)
  }

  if (funname == "MCP") { # total number of obs for MCP : n+3
    fit = mcp.fit(x, y, tau, lambda, a, weights = weights, warm = fit,
                  thresh = thresh, maxit = maxit, verbose = verbose)
  }

  fit = list(beta = fit, dim = dim(fit), lambda = lambda, df = sum(abs(fit) >0))
  class(fit) = "qcdfit"
  fit
}
