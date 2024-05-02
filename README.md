QCD : Coordinate Descent for Weighted Quantile RegressionDME
================

## Installing QCD package

``` r
install.packages("devtools")
library(devtools)

# check if this works
install_github("sangheekim96/QCD/QCD") # check if this works

# if it does not, we need to create a token due to 'private' repository
usethis::use_git_config(user.name = "someone", user.email = "someone@email.com")
usethis::create_github_token() # once the token is created, copy the Personal Access Token
credentials::set_github_pat() # run this one time and paste PAT
install_github("sangheekim96/QCD/QCD") # then try this

library(QCD)
```

## Description of `QCD` package

The `QCD` package consists of 7 functions.

- generate.data
- lasso.fit
- scad.fit
- mcp.fit
- qcd.fit
- qcd.path
- rmse

Each function has a help page!
