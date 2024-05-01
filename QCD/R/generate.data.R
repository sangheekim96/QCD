## FUNCTION : generate Peng data

library(mvtnorm)
generate.data = function(n, p, signal = 1) {
  indx = c(6, 12, 15, 20)
  true_beta = rep(0, 150)
  true_beta[indx] = signal
  true_beta[1] = 0.7*qnorm(0.3, mean=0, sd=1)
  Sigma = 0.5^abs(outer(1:p,1:p,'-'))
  X = rmvnorm(n, mean=rep(0,p), sigma=Sigma)
  epsilon = rnorm(n)
  Y = X[,6]+X[,12]+X[,15]+X[,20]+0.7*pnorm(X[,1])*epsilon

  return(list(Y=Y, X=X, true_beta=true_beta))
}


