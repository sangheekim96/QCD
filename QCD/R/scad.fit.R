#' Solve SCAD penalized quantile regression by coordinate descent for a single value of lambda
#'
#' @param x n x p design matrix X
#' @param y n x 1 vector of response variables Y
#' @param tau quantile
#' @param lambda single lambda value
#' @param a threshold parameter that adjusts constant penalty part (a > 2 for SCAD)
#' @param weights can input different weight for each p coefficients (default is no weights)
#' @param warm whether warm-start will be used (default is NULL). do not use this argument when only running scad.fit(). it is used when constructing regularization path with qcd.path() later.
#' @param thresh threshold of checking whether the coefficients converged
#' @param maxit maximum iteration for convergence
#' @param verbose  whether the iteration number will be printed (verbose = TRUE will print the iteration)
#'
#' @return beta  A n x 1 matrix of coefficients, stored in sparse matrix format
#' @return dim  dimension of coefficient matrix
#' @return lambda  lambda value used
#' @return df  number of nonzero coefficients
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
#' qr.scad = scad.fit(x = x, y = y, tau = 0.5, lambda = 0.8, a = 2.2,
#' weights = NULL, warm = NULL, thresh = 1e-06,
#' maxit = 10000, verbose = TRUE)
scad.fit <- function(x, y, tau, lambda, a, weights = NULL, warm = NULL,
                    thresh = thresh, maxit = maxit, verbose = FALSE) {

  n = nrow(x)
  p = ncol(x)

  if (is.null(weights))  weights = c(rep(1,p))
  else weights = weights

  if (is.null(warm))  beta = c(rep(0,p)) # initialize betas
  else beta = warm  # use previous betahat as a warm start for next lambda

  betahat = matrix(rep(NA, p*maxit), ncol = p) # matrix to store the final betahats
  r = matrix(0, nrow = n+5, ncol = p)  # total number of obs for SCAD : n+5

  if (a < 2) stop("a should be larger than 2") # SCAD requires a > 2

  for (l in 1:maxit) {

    if (verbose == TRUE) {
      print(l)
      print(beta)
    }

    for (j in 1:p) {
      r[1:n,j] = (y - x[,-j] %*% matrix(beta[-j], ncol = 1)) / x[,j]
      r[(n+1):(n+5),j] = c(0,-a*lambda, -lambda, lambda, a*lambda)
      w = c(weights, rep(0,5))
      x_added = c(x[,j], rep(0,5))
      taus = c(ifelse(x_added[1:n]>0, tau, 1-tau), rep(0,5))

      or = order(r[,j])
      newr = r[,j][or]
      newx = x_added[or]
      newtaus = taus[or]
      neww = w[or]

      Si_prime = neww * abs(newx) * newtaus
      wx = neww * abs(newx) / n
      S_prime = -sum(Si_prime) / n
      indx0 = max(which(newr == 0)) # find r=0
      indx_al1 = max(which(newr == -a*lambda)) # find r=-a*lambda
      indx_l1 = max(which(newr == -lambda)) # find r=-lambda
      indx_al2 = max(which(newr == a*lambda)) # find r=a*lambda
      indx_l2 = max(which(newr == lambda)) # find r=lambda

      base = cumsum(c(S_prime, wx)) # no boost for all points
      boost1 = (-a*lambda - newr[indx_al1:(indx_l1-1)]) / (a-1) # r between -a*lambda and -lambda gets line boost
      boost2 = rep(-lambda, indx0-indx_l1) # r between -lambda and 0 gets -lambda boost
      boost3 = rep(lambda, indx_l2 - indx0) # r between 0 and lambda gets lambda boost
      boost4 = (a*lambda - newr[(indx_l2+1):indx_al2]) / (a-1) # r between lambda and a*lambda gets line boost
      gradient = c(base[1:indx_al1],
                   base[(indx_al1+1):indx_al2] + c(boost1, boost2, boost3, boost4),
                   base[(indx_al2):length(base)])
      indx = min(which(gradient >= 0))
      betahat[l,j] = newr[indx-1]
      beta[j]=betahat[l,j]
    }

    cur_norm = norm(betahat[l,], type="2")
    old_norm = norm(betahat[l-1,], type="2")

    if (abs(cur_norm - old_norm) < thresh) { # stopping rule
      break
    }
  }

  out = as.numeric(tail(na.omit(betahat), 1)) # store the last betahat that converged

  return(out)
}
