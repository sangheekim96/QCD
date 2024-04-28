# QCD : Coordinate Descent for Weighted Quantile Regression

- glmnet.path() uses glmnet.fit() to optimize for single lambda, inside glmnet.fit() elnet.fit() coded in C++ is used to ensure fast computation
- For us, qcd.path() uses qcd.fit() to optimize for single lambda, and lasso.fit() / scad.fit() / mcp.fit() is used inside qcd.fit()
- We need to transform lasso.fit() / scad.fit() / mcp.fit() into Fortran

  *** Note : when running simulation, it takes longer when I use the above structure than when I use the original warm-start 

1. qcd.path() : input grid of lambda's, construct regularization path, use warm-start (deafult)
2. qcd.fit() : solve regularized QR for single lambda, we can select which penalty to use (LASSO, SCAD, MCP)
   - lasso.fit() : use lasso penalty
   - scad.fit() : use scad penalty
   - mcp.fit() : use mcp penalty

