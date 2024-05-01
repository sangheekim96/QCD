lasso.fit = function(x, y, tau, lambda, weights = NULL, warm = NULL,
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
      wx = neww * abs(newx) / n
      S_prime = -sum(Si_prime) / n
      indx0 = max(which(newr == 0)) # find the index when r=0, if multiple, select the last one

      base = cumsum(c(S_prime, wx))
      boost = c(rep(-lambda,indx0), rep(lambda, length(base)-indx0))
      gradient = base + boost
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
