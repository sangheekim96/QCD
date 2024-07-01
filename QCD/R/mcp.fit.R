#' Solve MCP penalized quantile regression by coordinate descent for a single value of lambda
#'
#' @param x n x p design matrix X
#' @param y n x 1 vector of response variables Y
#' @param tau quantile
#' @param lambda single lambda value
#' @param a threshold parameter that adjusts constant penalty part (a > 1 for MCP)
#' @param weights can input different weight for each p coefficients (default is no weights)
#' @param warm whether warm-start will be used (default is NULL). do not use this argument when only running mcp.fit(). it is used when constructing regularization path with qcd.path() later.
#' @param thresh threshold of checking whether the coefficients converged
#' @param maxit maximum iteration for convergence
#' @param verbose  whether the iteration number will be printed (verbose = TRUE will print the iteration)
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
#' qr.mcp = mcp.fit(x = x, y = y, tau = 0.5, lambda = 0.8, a = 2.2,
#' weights = NULL, warm = NULL, thresh = 1e-06,
#' maxit = 10000, verbose = TRUE)
mcp.fit <- function(x, y, tau, lambda, a, weights = NULL, warm = NULL,
                   thresh = thresh, maxit = maxit, verbose = FALSE) {

  n = nrow(x)
  p = ncol(x)

  if (is.null(weights))  weights = c(rep(1,p))
  else weights = weights

  if (is.null(warm))  beta = c(rep(0,p)) # initialize betas
  else beta = warm  # use previous betahat as a warm start for next lambda

  betahat = matrix(rep(NA, p*maxit), ncol = p) # matrix to store the final betahats
  r = matrix(0, nrow = n+3, ncol = p)  # total number of obs for MCP : n+3

  if (a < 1) stop("a should be larger than 1") # MCP requires a > 1

  for (l in 1:maxit) {
    if (verbose == TRUE) {
      print(l)
      print(beta)
    }

    for (j in 1:p) {
      r[1:n,j] = (y - x[,-j] %*% matrix(beta[-j], ncol = 1))/x[,j]
      r[(n+1):(n+3),j] = c(0, -a*lambda, a*lambda)
      w = c(weights, rep(0,3))
      x_added = c(x[,j], rep(0,3))
      taus = c(ifelse(x_added[1:n]>0, tau, 1-tau), rep(0,3))

      or = order(r[,j])
      newr = r[,j][or]
      newx = x_added[or]
      newtaus = taus[or]
      neww = w[or]

      Si_prime = neww * abs(newx) * newtaus
      wx = neww * abs(newx)
      S_prime = -sum(Si_prime)
      indx0 = max(which(newr == 0)) # find r=0
      indx_al1 = max(which(newr == -a*lambda)) # find r=-a*lambda
      indx_al2 = max(which(newr == a*lambda)) # find r=a*lambda

      base = cumsum(c(S_prime, wx))
      boost1 = -lambda - (newr[indx_al1:(indx0-1)] / a) # r between -a*lambda and 0 gets -lambda-beta/a boost
      boost2 = lambda - (newr[indx0:(indx_al2-1)] / a) # r between 0 and a*lambda gets lambda-beta/a boost
      gradient = c(base[1:indx_al1],
                   base[(indx_al1+1):indx_al2] + c(boost1, boost2),
                   base[(indx_al2):length(base)])
      indx = min(which(gradient >= 0))
      betahat[l,j] = newr[indx-1]
      beta[j]=betahat[l,j]
    }

    cur_norm = norm(betahat[l,], type="2")
    old_norm = norm(betahat[l-1,], type="2")

    if (abs(cur_norm - old_norm) < thresh) {
      break
    }
  }

  out = as.numeric(tail(na.omit(betahat), 1)) # store the last betahat that converged

  return(out)
}
