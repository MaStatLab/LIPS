Scalable Bayesian model averaging through local information propagation
================================

This package the LIPS algorithm for scalable Bayesian model averaging on linear regression models. This package uses the C++ library smctc for implementing the sequantial Monte Carlo samplers. 

### Install
The package can be installed on Linux and Mac using `devtools`:

```S
install.packages('devtools')
library('devtools')
devtools::install_github('MaStatLab/LIPS')
```

### Use
There are three functions in this package, and their descriptions are provided in the help files

```S
ans = lips(Y,X)
ans = lips_hd(Y,X)
ans = lips_hd_pred(Y,X)
```

### Reference
Ma L. and Soriano J. (2015). Scalable Bayesian model averaging through local information propagation. Journal of the American Statistical Association, Theory and Methods. Vol. 110, No. 510, 795-809.
