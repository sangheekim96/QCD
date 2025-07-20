QCD : Pathwise Coordinate Descent for High dimensional Penalized Quantile Regression
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

The purpose of this section is to give users a general sense of the package regarding QCD. We will briefly go over the main functions, basic operations and outputs. `QCD` can be loaded using the `library` command:

```{r}
library(QCD)
```

We create a data set for illustration:

```{r}
set.seed(1)
data <- generate.data(n = 30, p = 30, tau = 0.5, signal = 1, autocorrelation = 0.5)
x <- data$X
y <- data$Y
```

This function generates simulation data based on Peng and Wang (2015). Note that in our work, we do not consider intercept and we can change the value of true beta through signal argument. 

With single value of $\lambda$, we can fit the $\ell_1$ penalized quantile regression using 'qcd.lasso.fit'.

```{r}
qr.lasso <- qcd.lasso.fit(x = x, y = y, 
                          tau = 0.5, lambda = 0.01,
                          thresh = 1e-06, maxit = 100)
qr.lasso
```

We can fit the SCAD and MCP penalized quantile regression using 'qcd.scad.fit' and 'qcd.mcp.fit'. These functions are experimental.

```{r}
qr.scad <- qcd.scad.fit(x = x, y = y, 
                        tau = 0.5, lambda = 0.01,
                        a = 2.2,
                        thresh = 1e-06, maxit = 100)
qr.scad
```

```{r}
qr.mcp <- qcd.mcp.fit(x = x, y = y, 
                      tau = 0.5, lambda = 0.01,
                      a = 2.2,
                      thresh = 1e-06, maxit = 100)
qr.mcp
```

Users are encouraged to utilize 'qcd.path' function to solve penalized quantile regression problems. A pathwise scheme 'warm start' is used as a default.

```{r}
## Create lambda grid
upper <- 2; lower <- -6
lambda <-  2^seq(upper, lower, by = -0.2)

## warm start version
qr.lasso.warm = qcd.path(x = x, y = y, tau = 0.5,
                         funname = "LASSO", lambda = lambda,
                         nudge = FALSE, 
                         thresh = 1e-06, maxit = 10000)

qr.lasso.warm

## warm start and nudge version
set.seed(1)
qr.lasso.warm.nudge = qcd.path(x = x, y = y, tau = 0.5,
                               funname = "LASSO", lambda = lambda,
                               nudge = TRUE, nudgesd = 0.2, 
                               thresh = 1e-06, maxit = 10000)

qr.lasso.warm.nudge
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
