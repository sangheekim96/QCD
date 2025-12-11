QCD : A Pathwise Coordinate Descent Algorithm for LASSO Penalized Quantile Regression
================

## Introduction

We provide a `QCD` package that solves penalized quantile regression via exact coordinate descent method. The penalties considered are LASSO, SCAD, and MCP. Note that QCD algorithms for SCAD and MCP are experimental. Details may be found in Kim and Basu ([2025](#ref-qcd)) and in the vignette. Please email Sanghee Kim <sk2689@cornell.edu> if any bugs/errors have been discovered.

This file describes basic usage of functions related to $\ell_1$ penalized quantile regression model in `QCD` package in R.

`QCD` mainly solves the following problem. Given data points $(x_1, y_1), \ldots, (x_n, y_n)$, where $y_i \in \mathbb{R}$ is a numerical response variable, and $x_i \in \mathbb{R}^p$ is a $p$-dimensional covariate,

<p>
  $$
  \arg \min_{\beta \in \mathbb{R}^p} \sum_{i=1}^n \rho_\tau \left(y_i - x_i^\top \beta \right)
  + \lambda \sum_{j=1}^p \left| \beta_j \right|, \quad
  $$
</p> 
where &rho;<sub>&tau;</sub>(u) := u (&tau; - 𝐈(u < 0)) is the check loss function, and &lambda; is a penalty parameter to be chosen in a data-driven fashion. &lambda; could be a grid of values covering the entire range of possible solutions.

## Installation

You can download the `QCD` package from Github.

```{r, warning=FALSE, results='hide',message=FALSE}
library(devtools)
install_github("sangheekim96/QCD")
```

## Quick Start: QCD
The purpose of this section is to give users a general sense of the package regarding QCD. We will briefly go over the main functions, basic operations and outputs. `QCD` can be loaded using the `library` command in the standard way:


```{r}
library(QCD)
```
QCD solves for the exact solution of penalized quantile regression under high-dimensional settings.

### High-Dimensional Data Generation
The function generate.data() simulates data for penalized quantile regression, based on Peng and Wang (2015). It creates a sparse linear model with correlated covariates and heteroskedastic noise. User can freely adjust signal-to-noise ratio and autocorrelation rate in the covariance function. Note that in our work, we do not consider intercept and we can change the value of true beta through signal argument.

We generate a high-dimensional dataset for illustration, where the number of predictors exceeds the number of observations with n = 50, p = 150. Reader can also specify different configurations of data generating process in the following code chunk.
```{r}
n <- 50
p <- 150
set.seed(9)
data <- generate.data(n = n, p = p, tau = 0.3,
                       signal = 1, autocorrelation = 0.5)
x <- data$X
y <- data$Y
```


### QCD Under Lasso Penalty
With single value of $\lambda$, we can fit the $\ell_1$ penalized quantile regression using 'qcd.lasso.fit'.
```{r, results=FALSE}
qr.lasso <- qcd.lasso.fit(x = x, y = y,
tau = 0.3, lambda = 0.01,
thresh = 1e-06, maxit = 100)
```
This function uses exact coordinate-wise minimization for the quantile check loss, which attains similar accuracy over approximate methods. It is particularly suitable for high-dimensional problems (i.e., when `p >> n`), where traditional solvers like linear programming (LP) become computationally expensive.

### QCD Under Non-Convex Penalties

We can fit the SCAD and MCP penalized quantile regression using ‘qcd.scad.fit’ and ‘qcd.mcp.fit’ with single value of $\lambda$. These functions are experimental.

```{r, results=FALSE}
qr.scad <- qcd.scad.fit(x = x, y = y, 
                        tau = 0.3, lambda = 0.01,
                        a = 2.2,
                        thresh = 1e-06, maxit = 100)
qr.scad
```

```{r, results=FALSE}
qr.mcp <- qcd.mcp.fit(x = x, y = y, 
                      tau = 0.3, lambda = 0.01,
                      a = 2.2,
                      thresh = 1e-06, maxit = 100)
qr.mcp
```

### Across a grid of lambdas
The `qcd.path()` function solves penalized quantile regression over a grid of `lambda` values using **pathwise coordinate descent algorithm**. It supports LASSO, SCAD, and MCP penalties, and can optionally apply warm start and nudge to improve numerical stability and computation time along the solution path.

```{r, results=FALSE}
## Create lambda grid
upper <- 5; lower <- -7
lambda.grid <-  2^seq(upper, lower, by = -0.2)

## warm start version
qr.lasso.warm = qcd.path(x = x, y = y, tau = 0.3,
                         funname = "LASSO", lambda = lambda.grid,
                         nudge = FALSE, 
                         thresh = 1e-06, maxit = 1000)

qr.lasso.warm

## warm start and nudge version
set.seed(1)
qr.lasso.warm.nudge = qcd.path(x = x, y = y, tau = 0.3,
                               funname = "LASSO", lambda = lambda.grid,
                               nudge = TRUE, nudgesd = 0.01, 
                               thresh = 1e-06, maxit = 1000)

qr.lasso.warm.nudge
```

### Comparison of RMSE path
We will evaluate the estimation accuracy of QCD algorithm using three different strategies: vanilla, warm start, and warm nudge. We will measure the performance using RMSE.

```{r}
# true coefficient
true_beta <- data$true_beta

# Helper function for Relative Squared Error
calc_rmse <- function(beta_hat, beta_true) {
  sum((beta_hat - beta_true)^2) / sum(beta_true^2)
}

# Vanilla
rmse.vanilla <- numeric(length(lambda.grid))
b0 <- rep(0, p) 

for (l in seq_along(lambda.grid)) {
  fit <- qcd.lasso.fit(x = x, y = y, tau = 0.3, 
                       lambda = lambda.grid[l], 
                       warm = NULL, # No warm start
                       thresh = 1e-4, maxit = 100)
  
  rmse.vanilla[l] <- calc_rmse(fit$beta, true_beta)
}

# Warm Start
fit.warm <- qcd.path(x = x, y = y, tau = 0.3,
                     funname = "LASSO",
                     lambda = lambda.grid,
                     nudge = FALSE,
                     thresh = 1e-4, maxit = 100)

rmse.warm <- apply(fit.warm$beta, 2, calc_rmse, beta_true = true_beta)

# Warm Nudge
set.seed(1)
fit.nudge <- qcd.path(x = x, y = y, tau = 0.3,
                      funname = "LASSO",
                      lambda = lambda.grid,
                      nudge = TRUE, nudgesd = 0.01,
                      thresh = 1e-4, maxit = 100)

rmse.warm.nudge <- apply(fit.nudge$beta, 2, calc_rmse, beta_true = true_beta)

# Plot the RMSE path
QCD_df <- data.frame(
  lambda = log2(lambda.grid),
  rmse = c(rmse.vanilla, rmse.warm, rmse.warm.nudge),
  Method = rep(c('Vanilla', 'Warm', 'Warm nudge'), each = length(lambda.grid))
)

library(ggplot2)
ggplot(QCD_df, aes(x = lambda, y = rmse)) +
  geom_line(aes(colour = Method), linewidth = 1) +
  labs(x = expression(log[2](lambda)), y = "RMSE") +
  theme_classic() +
  scale_colour_manual(values = c("black", "blue", "red")) + 
  coord_cartesian(xlim = c(-7, 5), ylim = c(0, 1)) 
```
From above RMSE plot, vanilla QCD displays bumpy regularization path. Once warm start is implemented, the regularization path smoothens but the RMSE remains constant after a specific $\lambda$, indicating it is stuck at a suboptimal point. If a nudge is added to the warm start, the regularization path shows a U shape, suggesting that nudge helps our QCD algorithm converge to an optimal point.



### Runtime Comparison with LP and QICD
We compare the runtime of QCD with the benchmark methods: LP and Quantile Iterative Coordinate Descent (QICD). QCD runs faster than LP and QICD.

```{r, message=FALSE, warning=FALSE, eval = requireNamespace("quantreg", quietly = TRUE)}
library(quantreg)
library(QICD)

p <- 300
n <- 50
tau <- 0.3
upper <- 5; lower <- -5
lambda.grid <- 2^seq(upper, lower, by = -0.2)

data <- generate.data(p = p, n = n, signal = 1, tau = tau)
x <- data$X
y <- data$Y
  
# Method 1: QICD 
t.qicd <- NA
ptm <- proc.time()
for (l in seq_along(lambda.grid)) {
  fit.qicd <- QICD(x = x, y = y, tau = tau, 
                   lambda = lambda.grid[l], intercept = FALSE, funname = 'lasso')
}
t.qicd <- proc.time() - ptm
  
# Method 2: QCD
set.seed(9)
ptm <- proc.time()
fit.qcd <- qcd.path(x = x, y = y, tau = tau, 
                    funname = "LASSO", 
                    lambda = lambda.grid,
                    nudge = TRUE, nudgesd = 0.01, # Warm start with Nudge
                    thresh = 1e-4, maxit = 100)
t.qcd <- proc.time() - ptm

# Method 3: Linear Programming (quantreg)
ptm <- proc.time()
for (l in seq_along(lambda.grid)) {
  fit.lp <- rq.fit.lasso(x = x, y = y, tau = tau, 
                         lambda = rep(lambda.grid[l]*2, p))
}
t.lp <- proc.time() - ptm
  
print(t.lp)
print(t.qicd)
print(t.qcd)
```

## References

<div id="refs" class="references">

<div id="ref-qcd">

Sanghee Kim and Sumanta Basu. 2025. 
“A Pathwise Coordinate Descent Algorithm for LASSO Penalized Quantile Regression.” *arXiv preprint arXiv:2502.12363*. 
<https://doi.org/10.48550/arXiv.2502.12363>.

<div id="refs" class="references">

<div id="ref-qicd">

Peng Bo and Lan Wang. 2015. 
“An iterative coordinate descent algorithm for high-dimensional nonconvex penalized quantile regression.” *Journal of Computational and Graphical Statistics* 24.3 (2015): 676-694. 
<https://doi.org/10.1080/10618600.2014.913516>.

</div>

<div id="refs" class="references">

<div id="ref-glmnet">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.

</div>
