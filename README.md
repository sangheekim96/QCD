QCD : Pathwise Coordinate Descent for High dimensional Penalized Quantile Regression
================

# QCD

We provide a QCD algorithm that solves penalized quantile regression problem through exact pathwise coordinate descent. The penalties considered are LASSO, SCAD, and MCP. Note that QCD algorithms for SCAD and MCP are experimental. Details may be found in Kim and Basu ([2025](#ref-qcd)) and in the vignette. Please email Sanghee Kim <sk2689@cornell.edu> if any bugs/errors have been discovered.

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
