#' Solve Lasso penalized quantile regression by coordinate descent for a single value of lambda
#'
#' @param x n x p design matrix X
#' @param y n x 1 vector of response variables Y
#' @param tau quantile
#' @param lambda single lambda value
#' @param weights can input different weight for each p coefficients (default is no weights)
#' @param warm whether warm-start will be used (default is NULL). do not use this argument when only running lasso.fit(). it is used when constructing regularization path with qcd.path() later.
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
#' qr.lasso = lasso.fit(x = x, y = y, tau = 0.5, lambda = 0.8,
#' weights = NULL, warm = NULL, thresh = 1e-06,
#' maxit = 10000, verbose = TRUE)

lasso.fit <- function(x, y, tau, lambda, weights = NULL, warm = NULL,
                     thresh = thresh, maxit = maxit, verbose = FALSE) {

  n = nrow(x)
  p = ncol(x)

  if (is.null(weights))  weights = c(rep(1,p))
  else weights = weights

  if (is.null(warm))  beta = c(rep(0,p)) # initialize betas
  else beta = warm  # use previous betahat as a warm start for next lambda

  betahat = matrix(rep(NA, p*maxit), ncol = p) # matrix to store the final betahats
  r = matrix(0, nrow = n+1, ncol = p)  # total number of obs for Lasso : n+1

  for (l in 1:maxit) {

    if (verbose == TRUE) {
      print(l)
    }

    for (j in 1:p) {  #find the gradient for each beta_j
      r[1:n,j] = (y - x[,-j] %*% matrix(beta[-j], ncol = 1)) / x[,j]
      r[n+1,j] = 0 # add r_n+1 = 0
      w = c(weights, 0)
      x_added = c(x[,j], 0)
      taus = c(ifelse(x_added[1:n]>0, tau, 1-tau), 0)

      or = order(r[,j]) # order the residuals to compute gradient within each range
      newr = r[,j][or]
      newx = x_added[or]
      newtaus = taus[or]
      neww = w[or]

      Si_prime = neww * abs(newx) * newtaus
      wx = neww * abs(newx)
      S_prime = -sum(Si_prime)
      indx0s = which(newr == 0)
      indx0 = max(indx0s) # find the index when r=0, if multiple, select the last one
      
      # KKT condition
      if ( abs(S_prime + sum(wx[1:(min(indx0s)-1)])) <= lambda ) {
        betahat[l,j] = 0
        beta[j] = 0
      }

      else {
        base = cumsum(c(S_prime, wx))
        boost = c(rep(-lambda,indx0), rep(lambda, length(base)-indx0))
        gradient = base + boost
        indx = min(which(gradient >= 0)) 
        betahat[l,j] = newr[indx-1] 
        beta[j] = betahat[l,j]
      }
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
