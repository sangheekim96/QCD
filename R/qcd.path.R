qcd.path <- function(x, y, tau, funname = c("LASSO", "SCAD", "MCP"), a, weights = NULL,
                     lambda = NULL, nlambda = 30, nudge = FALSE, nudgesd = 0.01, standardize = TRUE,
                     thresh = 1e-6, maxit = 100000, verbose = TRUE) {

  # store dimension
  np = dim(x)
  if (is.null(np) || (np[2] <= 1)) stop("x should be a matrix with 2 or more columns")
  n = as.integer(np[1]); p = as.integer(np[2])

  # get feature variable names
  vnames=colnames(x)
  if (is.null(vnames)) vnames = paste("V",seq(p),sep="")

  # check weights
  if (is.null(weights)) weights = rep(1,n)
  else if (length(weights) != n)
    stop(paste("Number of elements in weights (",length(weights),
               ") not equal to the number of rows of x (",n,")",sep=""))
  weights <- as.double(weights)

  # check standardize
  if (standardize) xs = scale(x)
  else xs = x

  # work out lambda values
  nlam = as.integer(nlambda)
  user_lambda = FALSE   # did user provide their own lambda values?
  if (is.null(lambda)) stop("lambda needs to be specified")
  else {  # user provided lambda values
    user_lambda = TRUE
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    ulam = as.double(sort(lambda, decreasing = TRUE))
    nlam = as.integer(length(lambda))
  }

  # check if warm start will be used and create the path
  path = matrix(NA, nrow=length(lambda), ncol=p) # create matrix to store the path for different lambdas

  fit = NULL
  #set.seed(9)
  for (k in 1:nlam) {
    cur_lambda = ulam[k]
    cat("Fitting lambda index", k, ":", round(log2(cur_lambda), 2), fill = TRUE) # print which lambda we are at

    betahat = matrix(rep(NA, p*maxit), ncol=p)
    fit = qcd.fit(x = xs, y, tau, weights, cur_lambda, a, funname = funname, warm = fit,
                  thresh = thresh, maxit = maxit, verbose = FALSE)
    path[k,] = fit$beta

    # If nudge = TRUE, we can choose the sd of rnorm() by nudgesd, recommend to set.seed before qcd.path()
    if (nudge == TRUE) fit = path[k,] + rnorm(p, mean = 0, sd = nudgesd)
    else fit = path[k,]
  }

  # output
  out = list(beta = path, dim = dim(path), lambda = ulam, nobs = n, df = rowSums(abs(path) >0))
  class(out) = "qcd.path"

  return(out)
}
