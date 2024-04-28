# Use coordinate descent to solve regularized QR for a single value of lambda - MCP

mcp.fit = function(x, y, tau, lambda, a, weights = weights, warm = NULL,
                   thresh = thresh, maxit = maxit, verbose = FALSE) { 
  
  n = nrow(x)
  p = ncol(x)
  
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
      wx = neww * abs(newx) / n
      S_prime = -sum(Si_prime) / n
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
