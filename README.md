## Installing all packages needed

You can download the `QCD` package by devtools::install_github function. `QCD` package contains  functions that solve regularized weighted quantile regression problem through coordinate descent.
The penalties are Lasso, SCAD, and MCP. One can also run a pathwise coordinate descent through function that uses warm-start.

```{r setup, eval=FALSE}
install.packages("devtools")
library(devtools)
install_github("sangheekim96/QCD/QCD")

library(QCD)
```


## Description of `QCD` package

The `QCD` package consists of 7 functions.

## generate.data

This function generates simulation data based on Peng and Wang (2015). Note that in our work, we do not consider intercept and we can change the value of true beta through signal argument (default is 1 as in Peng and Wang (2015) paper). 

### Arguments
- n : number of observations
- p : number of variables
- signal : signal of true beta's

### Value
- Y : simulated n x 1 Y vector
- X : simulated n x p X matrix
- true_beta : simulated true beta's

```{r setup1, eval=FALSE}
set.seed(1)
data = generate.data(n = 150, p = 150, signal = 1)
```

## lasso.fit

This function solves Lasso penalized quantile regression by coordinate descent for a single value of lambda. It has the similar structure as elnet.fit().

### Arguments
- x : n x p design matrix X
- y : n x 1 vector of response variables Y
- tau : quantile
- lambda : single lambda value
- weights : can input different weight for each p coefficients (default is no weights)
- warm : whether warm-start will be used (default is NULL). we do not use this argument when only running lasso.fit(), but it is used when constructing regularization path with qcd.path() later. 
- thresh : threshold of checking whether the coefficients converged
- maxit : maximum iteration for convergence
- verbose : whether the iteration number will be printed. verbose = TRUE will print the iteration.

### Value
- beta : A n x 1 matrix of coefficients, stored in sparse matrix format
- dim : dimension of coefficient matrix
- lambda : lambda value used
- df : number of nonzero coefficients


```{r setup2, eval=FALSE}
qr.lasso = lasso.fit(x = data$X, y = data$Y, tau = 0.5, lambda = 0.8, 
                     weights = NULL, warm = NULL, thresh = 1e-06, 
                     maxit = 10000, verbose = TRUE)
```


## scad.fit / macp.fit

This function solves SCAD/MCP penalized quantile regression by coordinate descent for a single value of lambda. It has the similar structure as elnet.fit(). All arguments are the same as lasso.fit() and add argument `a`.

### Arguments
- x : n x p design matrix X
- y : n x 1 vector of response variables Y
- tau : quantile
- lambda : single lambda value
- a : threshold parameter that adjusts constant penalty part (a > 2 for SCAD, a > 1 for MCP)
- weights : can input different weight for each p coefficients (default is no weights)
- warm : whether warm-start will be used (default is NULL). we do not use this argument when only running lasso.fit(), but it is used when constructing regularization path with qcd.path() later. 
- thresh : threshold of checking whether the coefficients converged
- maxit : maximum iteration for convergence
- verbose : whether the iteration number will be printed. verbose = TRUE will print the iteration.

### Value
- beta : A n x 1 matrix of coefficients, stored in sparse matrix format
- dim : dimension of coefficient matrix
- lambda : lambda value used
- df : number of nonzero coefficients


```{r setup3, eval=FALSE}
qr.scad = scad.fit(x = data$X, y = data$Y, tau = 0.5, lambda = 0.8, a = 2.2, 
                   weights = NULL, warm = NULL, thresh = 1e-06, 
                   maxit = 10000, verbose = TRUE)

qr.mcp = mcp.fit(x = data$X, y = data$Y, tau = 0.5, lambda = 0.8, a = 2.2, 
                 weights = NULL, warm = NULL, thresh = 1e-06, 
                 maxit = 10000, verbose = TRUE)
```


## qcd.fit

This function is a general version where the user can choose which penalty they will use. It has the similar structure as glmnet.fit(). It implements lasso.fit(), scad.fit(), and mcp.fit() depending on user's choice. All arguments are the same as lasso.fit() and add argument `funname`.

### Arguments
- x : n x p design matrix X
- y : n x 1 vector of response variables Y
- tau : quantile
- lambda : single lambda value
- a : only needed when user chooses "SCAD" or "MCP". threshold parameter that adjusts constant penalty part (a > 2 for SCAD, a > 1 for MCP)
- funname : user can choose "LASSO", "SCAD", and "MCP" as the penalty
- weights : can input different weight for each p coefficients (default is no weights)
- warm : whether warm-start will be used (default is NULL). we do not use this argument when only running lasso.fit(), but it is used when constructing regularization path with qcd.path() later. 
- thresh : threshold of checking whether the coefficients converged
- maxit : maximum iteration for convergence
- verbose : whether the iteration number will be printed. verbose = TRUE will print the iteration.

### Value
- beta : A n x 1 matrix of coefficients, stored in sparse matrix format
- dim : dimension of coefficient matrix
- lambda : lambda value used
- df : number of nonzero coefficients


```{r setup4, eval=FALSE}
qr.lasso.choice = qcd.fit(x = data$X, y = data$Y, tau = 0.5, lambda = 0.8,
                          funname = "LASSO", weights = NULL, warm = NULL, 
                          thresh = 1e-06, maxit = 10000, verbose = TRUE)

qr.scad.choice = qcd.fit(x = data$X, y = data$Y, tau = 0.5, lambda = 0.8, a = 2.2,
                          funname = "SCAD", weights = NULL, warm = NULL, 
                          thresh = 1e-06, maxit = 10000, verbose = TRUE)

qr.mcp.choice = qcd.fit(x = data$X, y = data$Y, tau = 0.5, lambda = 0.8, a = 2.2,
                          funname = "MCP", weights = NULL, warm = NULL, 
                          thresh = 1e-06, maxit = 10000, verbose = TRUE)

```


## qcd.path

This function solves penalized quantile regression through pathwise coordinate descent. Warm-start is used as a default. qcd.path() uses qcd.fit() inside and denotes warm = fit to use warm-start. All arguments are the same as qcd.fit() except for `lambda` and `warm`.

### Arguments
- x : n x p design matrix X
- y : n x 1 vector of response variables Y
- tau : quantile
- funname : user can choose "LASSO", "SCAD", and "MCP" as the penalty
- a : only needed when user chooses "SCAD" or "MCP"
- weights : can input different weight for each p coefficients (default is no weights)
- lambda : vector of lambda gird should be input instead of a single value
- nlambda : number of lambda grid. it will be automatically calculated if there is no input.
- nudge : whether nudge will be added when changing lambda (default is nudge = FALSE). recommend set.seed() before the function if nudge is used.
- nudgesd : user can choose the standard deviation of the random nudge (default is 0.01)
- standardize : whether X will be scaled (default is standardize = TRUE)
- thresh : threshold of checking whether the coefficients converged
- maxit : maximum iteration for convergence
- verbose : whether the iteration number will be printed. verbose = TRUE will print the iteration.

### Value
- beta : A n x 1 matrix of coefficients, stored in sparse matrix format
- dim : dimension of coefficient matrix
- lambda : lambda value used
- nobs : number of observations
- df : number of nonzero coefficients for each lambda


```{r setup5, eval=FALSE}
# create lambda vector
upper = 2; lower = -6
lambda =  2^seq(upper, lower, by = -0.2)


# pathwise coordinate dewcent
qr.lasso.warm = qcd.path(x = data$X, y = data$Y, tau = 0.5, funname = "LASSO",
                         weights = NULL, lambda = lambda, nlambda = length(lambda),
                         nudge = FALSE, standardize = TRUE, 
                         thresh = 1e-06, maxit = 10000, verbose = TRUE)

set.seed(1)
qr.lasso.warm.nudge = qcd.path(x = data$X, y = data$Y, tau = 0.5, funname = "LASSO",
                               weights = NULL, lambda = lambda, 
                               nlambda = length(lambda), 
                               nudge = TRUE, nudgesd = 0.01, standardize = TRUE, 
                               thresh = 1e-06, maxit = 10000, verbose = TRUE)

qr.scad.warm = qcd.path(x = data$X, y = data$Y, tau = 0.5, funname = "SCAD", a = 2.2,
                               weights = NULL, lambda = lambda, 
                               nlambda = length(lambda), 
                               nudge = TRUE, nudgesd = 0.01, standardize = TRUE, 
                               thresh = 1e-06, maxit = 10000, verbose = TRUE)

```


## rmse

This function calculates RMSE of the estimated coefficients. The 

### Arguments
- beta : estimated beta's
- truebeta : true beta's
- lambda : grid of lambda's


```{r setup6, eval=FALSE}
lasso.warm.rmse = rmse(qr.lasso.warm$beta, data$true_beta, lambda)
```


