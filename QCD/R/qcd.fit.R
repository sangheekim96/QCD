# Use coordinate descent to solve regularized QR for a single value of lambda

qcd.fit = function(x, y, tau, weights, lambda, a, funname = c("LASSO", "SCAD", "MCP"),
                   warm = NULL, thresh = 1e-6, maxit = 1000000, verbose = FALSE) {

  if (is.null(warm)) fit = NULL
  else fit = warm # warm start input

  # run coordinate descent for different penalties
  if (funname == "LASSO") {
    fit = lasso.fit(x, y, tau, lambda, weights = weights, warm = fit,
                    thresh = thresh, maxit = maxit, verbose = verbose)
  }

  if (funname == "SCAD") { # total number of obs for SCAD : n+5
    fit = scad.fit(x, y, tau, lambda, a, weights = weights, warm = fit,
                   thresh = thresh, maxit = maxit, verbose = verbose)
  }

  if (funname == "MCP") { # total number of obs for MCP : n+3
    fit = mcp.fit(x, y, tau, lambda, a, weights = weights, warm = fit,
                  thresh = thresh, maxit = maxit, verbose = verbose)
  }

  fit = list(beta = fit, dim = dim(fit), lambda = lambda, df = sum(abs(fit) >0))
  class(fit) = "qcdfit"
  fit
}
