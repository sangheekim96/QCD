#' Solve MCP penalized quantile regression by coordinate descent for a single value of lambda
#'
#' @param x \code{n x p} design matrix X.
#' @param y \code{n x 1} vector of response variables Y.
#' @param tau Quantile value between 0 and 1.
#' @param lambda A single lambda value.
#' @param a Threshold parameter that adjusts constant penalty part. (a > 1 for MCP)
#' @param warm Whether warm start will be used.
#' Do not use this argument when only running \code{qcd.lasso.fit}.
#' It is used when constructing regularization path with \code{qcd.path} later.
#' @param thresh Threshold of checking whether the coefficients converged.
#' @param maxit Maximum iteration for convergence.
#'
#' @return \item{beta}{A \code{n x 1} matrix of coefficients.}
#' @return \item{dim}{Dimension of coefficient vector.}
#' @return \item{lambda}{Lambda value used.}
#' @return \item{df}{Number of nonzero coefficients.}
#' @import utils
#' @import stats
#' @export
#'
#' @examples
#' ## Create sample data set
#' n <- 30; p <- 30
#' x <- array(rnorm(n*p), c(n,p))
#' for (j in 1:p){x[,j] = x[,j]/(norm(as.matrix(x[,j]), type="f")/sqrt(n))}
#' e <- rnorm(n)
#' b <- c(1, -1, rep(0, p-2))
#' y <- x %*% b + e
#'
#' ## Use QCD algorithm to solve SCAD penalized quantile regression with one lambda
#' qr.mcp <- qcd.mcp.fit(x = x, y = y, tau = 0.5, lambda = 0.1, a = 2.2,
#' thresh = 1e-06, maxit = 100)
#'
#'
qcd.mcp.fit <- function(x, y, tau,
                         lambda,
                         a,
                         thresh = 1e-6,
                         maxit = 100000,
                         warm = NULL) {

  nobs <- nrow(x)
  nvars <- ncol(x)

  if (is.null(warm))  beta <- c(rep(0, nvars)) # initialize betas
  else beta <- warm  # use previous betahat as a warm start for next lambda

  betahat <- matrix(rep(NA, nvars*maxit), ncol = nvars) # matrix to store the final betahats
  r <- matrix(0, nrow = nobs+3, ncol = nvars)  # total number of obs for MCP : n+3

  if (a < 1) stop("a should be larger than 1") # MCP requires a > 1

  for (l in 1:maxit) {

    for (j in 1:nvars) {
      r[1:nobs,j] <- (y - x[,-j] %*% matrix(beta[-j], ncol = 1))/x[,j]
      r[(nobs+1):(nobs+3),j] <- c(0, -a*lambda, a*lambda)
      w <- c(rep(1, nvars), rep(0,3))
      x_added <- c(x[,j], rep(0,3))
      taus <- c(ifelse(x_added[1:nobs]>0, tau, 1-tau), rep(0,3))

      or <- order(r[,j])
      newr <- r[,j][or]
      newx <- x_added[or]
      newtaus <- taus[or]
      neww <- w[or]

      Si_prime <- neww * abs(newx) * newtaus
      wx <- neww * abs(newx) / nobs
      S_prime <- -sum(Si_prime) / nobs
      indx0 <- max(which(newr == 0)) # find r=0
      indx_al1 <- max(which(newr == -a*lambda)) # find r=-a*lambda
      indx_al2 <- max(which(newr == a*lambda)) # find r=a*lambda

      base <- cumsum(c(S_prime, wx))
      boost1 <- -lambda - (newr[indx_al1:(indx0-1)] / a) # r between -a*lambda and 0 gets -lambda-beta/a boost
      boost2 <- lambda - (newr[indx0:(indx_al2-1)] / a) # r between 0 and a*lambda gets lambda-beta/a boost
      gradient <- c(base[1:indx_al1],
                    base[(indx_al1+1):indx_al2] + c(boost1, boost2),
                    base[(indx_al2):length(base)])
      indx <- min(which(gradient >= 0))
      betahat[l,j] <- newr[indx-1]
      beta[j] <- betahat[l,j]
    }

    cur_norm <- norm(betahat[l,], type="2")
    old_norm <- norm(betahat[l-1,], type="2")

    if (abs(cur_norm - old_norm) < thresh) {
      break
    }
  }

  beta <- as.numeric(tail(na.omit(betahat), 1)) # store the last betahat that converged
  out <- list(
    beta = beta,
    df = sum(abs(beta) > 0),
    lambda = lambda)

  class(out) = "qcdmcpfit"
  return(out)
}
