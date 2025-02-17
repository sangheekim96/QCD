#' Solve LASSO penalized quantile regression by coordinate descent for a single value of lambda
#'
#' @param x \code{n x p} design matrix X.
#' @param y \code{n x 1} vector of response variables Y.
#' @param tau Quantile value between 0 and 1.
#' @param lambda A single lambda value.
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
#'
#' @useDynLib QCD, .registration = TRUE
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
#' ## Use QCD algorithm to solve LASSO penalized quantile regression with one lambda
#' qr.lasso <- qcd.lasso.fit(x = x, y = y, tau = 0.5, lambda = 0.1,
#' thresh = 1e-06, maxit = 100)
#'
#'
qcd.lasso.fit <- function(x, y, tau,
                          lambda,
                          thresh = 1e-6,
                          maxit = 100000,
                          warm = NULL) {

  nobs <- nrow(x)
  nvars <- ncol(x)

  if (is.null(warm))  beta <- c(rep(0, nvars)) # initialize betas
  else beta <- warm  # use previous betahat as a warm start for next lambda

  lasso_fit <- .Fortran("qcdwarm",
                        x = matrix(as.double(x), nrow = nobs, ncol = nvars),
                        y = as.double(y),
                        n = as.integer(nobs),
                        p = as.integer(nvars),
                        tau = as.double(tau),
                        lambda = as.double(lambda),
                        b0 = as.double(beta),
                        bhat = as.double(rep(0, nvars)),
                        it = as.integer(1),
                        mt = as.integer(maxit),
                        tl = as.double(thresh)
                        )

  beta <- lasso_fit$bhat
  out <- list(beta = beta,
              df = sum(abs(lasso_fit$bhat) > 0),
              lambda = lasso_fit$lambda)

  class(out) = "qcdlassofit"
  return(out)
}
