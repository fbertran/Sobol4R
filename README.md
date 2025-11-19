<!-- README.md is generated from README.Rmd. Please edit that file -->


# Sobol4R <img src="man/figures/logo.png" align="right" width="200"/>

# Sobol Indices for Models with Fixed and Stochastic Parameters
## Frédéric Bertrand, Elizaveta Logosha and Myriam Maumy-Bertrand

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/fbertran/Sobol4R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/Sobol4R/actions/workflows/R-CMD-check.yaml)
[![R-CMD-check](https://github.com/fbertran/Sobol4R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/Sobol4R/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/SobolR4)](https://cran.r-project.org/package=SobolR4)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/SobolR4)](https://cran.r-project.org/package=SobolR4)
[![GitHub Repo stars](https://img.shields.io/github/stars/fbertran/SobolR4?style=social)](https://github.com/fbertran/SobolR4/)
[![DOI](https://zenodo.org/badge/18437255.svg)](https://zenodo.org/badge/latestdoi/18437255)
[![Codecov test coverage](https://codecov.io/gh/fbertran/SobolR4/branch/master/graph/badge.svg)](https://codecov.io/gh/fbertran/SobolR4?branch=master)
<!-- badges: end -->

This site and the examples reproduced on it were created by F. Bertrand, E. Logosha and M. Maumy-Bertrand.

## Installation

You can install the latest version of the Sobol4R package from the [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("Sobol4R")
```

You can install the latest version of the Sobol4R package from [github](https://github.com) with:


```r
devtools::install_github("fbertran/Sobol4R")
```

## Exemples

### Context and non random case

#### Test case : the non-monotonic Sobol g-function

The method of sobol requires 2 samples (there are 2 variables, all following the uniform distribution on $[0,1$])


```r
library(Sobol4R)
set.seed(4669)
n <- 50000
#X1_1 <- data.frame(matrix(runif(8 * n), nrow = n))
data(X1_1, package="Sobol4R")
#X2_1 <- data.frame(matrix(runif(8 * n), nrow = n))
data(X2_1, package="Sobol4R")
```



```r
library(sensitivity)
if(requireNamespace("sensitivity", quietly = TRUE)){
sensitivity::sobol.fun
}
#> function (X) 
#> {
#>     a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
#>     y <- 1
#>     for (j in 1:8) {
#>         y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
#>     }
#>     y
#> }
#> <bytecode: 0x7f93346e7140>
#> <environment: namespace:sensitivity>
```


```r
library(sensitivity)
if(requireNamespace("sensitivity", quietly = TRUE)){
gensol1 <- sensitivity::sobol(model = NULL, X1 = X1_1, X2 = X2_1, order = 2, nboot = 100)
Y1 = sensitivity::sobol.fun(gensol1$X)
data(Y1, package="Sobol4R")
}
```



```r
if(requireNamespace("sensitivity", quietly = TRUE)){
x1 <- sensitivity::tell(gensol1, Y1)
print(x1)
ggplot2::ggplot(x1)
rm(gensol1,X1_1,X2_1)
}
#> 
#> Call:
#> sensitivity::sobol(model = NULL, X1 = X1_1, X2 = X2_1, order = 2,     nboot = 100)
#> 
#> Model runs: 1850000 
#> 
#> Sobol indices
#>            original          bias  std. error    min. c.i.  max. c.i.
#> X1     0.7174136542 -1.469403e-04 0.007167772  0.705217099 0.73144839
#> X2     0.1743523954 -6.880819e-04 0.009722451  0.154020771 0.19323394
#> X3     0.0172802083 -8.311862e-05 0.010790956 -0.005620065 0.03850099
#> X4     0.0007144226 -2.892190e-04 0.010765169 -0.022922480 0.02339221
#> X5    -0.0066421180 -1.810232e-04 0.010810575 -0.030803073 0.01620777
#> X6    -0.0066098421 -1.984479e-04 0.010804915 -0.030782553 0.01626296
#> X7    -0.0065132242 -1.879218e-04 0.010817631 -0.030633006 0.01622637
#> X8    -0.0066549163 -1.738606e-04 0.010812077 -0.030793484 0.01626333
#> X1*X2  0.0638877952  4.530138e-04 0.012620762  0.038237826 0.09116831
#> X1*X3  0.0146165648  2.944769e-04 0.011224011 -0.009379035 0.03945923
#> X1*X4  0.0090055130  2.863625e-04 0.010919358 -0.013772926 0.03387554
#> X1*X5  0.0066288190  1.955101e-04 0.010807333 -0.016301110 0.03068314
#> X1*X6  0.0066916868  1.829647e-04 0.010822398 -0.016218406 0.03082603
#> X1*X7  0.0066645835  1.739142e-04 0.010814017 -0.016218959 0.03084353
#> X1*X8  0.0065265988  1.886866e-04 0.010800668 -0.016252113 0.03062897
#> X2*X3  0.0082407818  2.384301e-04 0.010793489 -0.014622190 0.03230174
#> X2*X4  0.0071595658  1.859228e-04 0.010843721 -0.016295392 0.03162920
#> X2*X5  0.0067026359  1.868787e-04 0.010814049 -0.016164785 0.03086200
#> X2*X6  0.0067462441  1.855100e-04 0.010804006 -0.016066366 0.03082248
#> X2*X7  0.0066878725  1.817594e-04 0.010812661 -0.016207838 0.03078866
#> X2*X8  0.0066648267  1.929796e-04 0.010809582 -0.016219286 0.03074448
#> X3*X4  0.0067687716  2.156748e-04 0.010799560 -0.015997302 0.03067781
#> X3*X5  0.0067072641  1.847805e-04 0.010809585 -0.016107244 0.03080750
#> X3*X6  0.0066888007  1.880596e-04 0.010811562 -0.016168034 0.03079173
#> X3*X7  0.0066850717  1.840728e-04 0.010810248 -0.016166273 0.03078838
#> X3*X8  0.0066828702  1.885483e-04 0.010809713 -0.016155296 0.03076988
#> X4*X5  0.0066761722  1.876135e-04 0.010809274 -0.016181971 0.03077676
#> X4*X6  0.0066969601  1.864462e-04 0.010810141 -0.016152524 0.03079303
#> X4*X7  0.0066963855  1.885056e-04 0.010808795 -0.016143164 0.03079264
#> X4*X8  0.0066836048  1.858056e-04 0.010809312 -0.016169452 0.03078736
#> X5*X6  0.0066872687  1.867842e-04 0.010809534 -0.016158725 0.03078624
#> X5*X7  0.0066868618  1.866769e-04 0.010809430 -0.016156747 0.03078493
#> X5*X8  0.0066861003  1.867500e-04 0.010809562 -0.016160417 0.03078676
#> X6*X7  0.0066870754  1.866652e-04 0.010809440 -0.016157834 0.03078701
#> X6*X8  0.0066871628  1.868050e-04 0.010809377 -0.016159054 0.03078519
#> X7*X8  0.0066872274  1.866333e-04 0.010809436 -0.016157482 0.03078613
```


