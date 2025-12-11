#' solve penalized quantile regression by pathwise coordinate descent for a grid of lambda's with user chosen penalty
#'
#' @param x \code{n x p} design matrix X.
#' @param y \code{n x 1} vector of response variables Y.
#' @param tau Quantile value between 0 and 1.
#' @param funname User chosen "LASSO", "SCAD", and "MCP" as the penalty.
#' @param a Threshold parameter that adjusts constant penalty part for SCAD and MCP. (a > 2 for SCAD, a > 1 for MCP)
#' @param lambda Vector of lambda gird.
#' @param nudge Whether nudge will be added when changing lambda. (default is nudge = FALSE) Recommended to set.seed() before the function if nudge is used.
#' @param nudgesd User chosen standard deviation of the random nudge. (default is 0.01)
#' @param standardize Whether X will be scaled.
#' @param thresh Threshold of checking whether the coefficients converged.
#' @param maxit Maximum iteration for convergence.
#'
#' @return \item{beta}{A \code{p x length(lambda)} matrix of coefficients, stored in
#' sparse matrix format.}
#' @return \item{dim}{Dimension of coefficient matrix.}
#' @return \item{lambda}{The actual sequence of lambda values used.}
#' @return \item{nobs}{Number of observations.}
#' @return \item{df}{The number of nonzero coefficients for each value of lambda.}
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
#' ## Create lambda grid
#' upper <- 2; lower <- -5
#' lambda <-  2^seq(upper, lower, by = -0.2)
#'
#' ## Use QCD algorithm to solve penalized quantile regression with given lambda grid
#' qr.lasso.warm = qcd.path(x = x, y = y, tau = 0.5,
#' funname = "LASSO", lambda = lambda,
#' nudge = FALSE, thresh = 1e-06, maxit = 10000)
#'
#' set.seed(1)
#' qr.lasso.warm.nudge = qcd.path(x = x, y = y, tau = 0.5,
#' funname = "LASSO", lambda = lambda,
#' nudge = TRUE, nudgesd = 0.2, thresh = 1e-06, maxit = 10000)
#'
#'
qcd.path <- function(x, y, tau,
                     funname = c("LASSO", "SCAD", "MCP"),
                     a = 2.2,
                     lambda = NULL,
                     nudge = FALSE,
                     nudgesd = 0.01,
                     standardize = FALSE,
                     thresh = 1e-6,
                     maxit = 100000) {

  # store dimension
  np <- dim(x)
  if(is.null(np) || (np[2] <= 1)) {
    stop("x should be a matrix with 2 or more columns")
  }

  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])

  # get feature variable names
  vnames <- colnames(x)
  if(is.null(vnames)) {
    vnames <- paste("V",seq(nvars),sep="")
  }


  # standardize x if necessary
  if (standardize) {
    for (j in 1:nvars){
      x[,j] <- x[,j] / sqrt(mean(Mod(x[,j])^2))
    }
  }

  # check lambda values error
  user_lambda <- FALSE   # did user provide their own lambda values?

  if (is.null(lambda)) {
    stop("lambda values should be provided by the user")
  } else {  # user provided lambda values
    user_lambda <- TRUE
    if (!is.numeric(lambda))
      stop("lambda must be numeric")

    if (any(!is.finite(lambda)))
      stop("lambda must be finite")

    if (any(lambda < 0))
      stop("lambdas should be non-negative")

    ulam <- as.double(sort(lambda, decreasing = TRUE))
    nlam <- as.integer(length(lambda))

    # remove near-duplicates (relative tol ~ 1e-12)
    lamu <- lambda[ c(TRUE, diff(lambda) < -1e-12 * lambda[-length(lambda)]) ]
    if (length(lamu) < length(ulam)) {
      warning("Removed near-duplicate lambda values (within 1e-12 relative).")
    }
  }

  # for scad, mcp penalties, check if a is present
  if (funname %in% c("SCAD","MCP")) {
    if (is.null(a)) stop("'a' must be provided for SCAD/MCP.")
    if (!is.numeric(a) || length(a) != 1L || !is.finite(a))
      stop("'a' must be a single finite numeric value.")

    if (funname == "SCAD") {
      if (a <= 2) stop("SCAD: 'a' must be > 2")

    } else { # MCP case
      if (a <= 1) stop("MCP: 'a' (gamma) must be > 1.")
    }
  }


  beta <- matrix(0, nrow = nvars, ncol = nlam)
  fit <- NULL

  for (k in 1:nlam) {
    cur_lambda <- ulam[k]
    #cat(log2(cur_lambda), " ")

    if (funname == "LASSO") {
      fit <- qcd.lasso.fit(x, y, tau,
                           cur_lambda,
                           thresh = thresh,
                           maxit = maxit,
                           warm = fit)
    }
    if (funname == "SCAD") {
      fit <- qcd.scad.fit(x, y, tau,
                          cur_lambda,
                          a = a,
                          thresh = thresh,
                          maxit = maxit,
                          warm = fit)
    }
    if (funname == "MCP") {
      fit <- qcd.mcp.fit(x, y, tau,
                         cur_lambda,
                         a = a,
                         thresh = thresh,
                         maxit = maxit,
                         warm = fit)
    }


    beta[, k] <- fit$beta

    # If nudge = TRUE, we can choose the sd of rnorm() by nudgesd
    # recommend to set.seed before qcd.path()
    if (nudge == TRUE) {
      fit <- beta[, k] + rnorm(nvars, mean = 0, sd = nudgesd)
    } else fit <- beta[, k]
  }

  # output
  stepnames <- paste0("s", 0:(length(ulam) - 1))
  dimnames(beta) <- list(vnames, stepnames)

  out <- list()
  out$beta <- beta
  out$df <- colSums(abs(beta) > 0)
  out$dim <- dim(beta)
  out$lambda <- ulam
  out$nobs <- nobs
  class(out) <- "qcdpath"

  return(out)
}
