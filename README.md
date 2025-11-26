<!-- README.md is generated from README.Rmd. Please edit that file -->


# Sobol4R, Sobol Indices for Models with Fixed and Stochastic Parameters <img src="man/figures/Sobol4R_hex.svg" align="right" width="200"/>

## Frédéric Bertrand, Elizaveta Logosha and Myriam Maumy-Bertrand

<!-- badges: start -->
[![R-CMD-check](https://github.com/fbertran/Sobol4R/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fbertran/Sobol4R/actions/workflows/R-CMD-check.yaml)
[![R-hub](https://github.com/fbertran/Sobol4R/actions/workflows/rhub.yaml/badge.svg)](https://github.com/fbertran/Sobol4R/actions/workflows/rhub.yaml)
<!-- badges: end -->

Tools to design experiments, compute Sobol sensitivity indices,
    and summarise stochastic responses inspired by the strategy described by
    Zhu et Sudret (2021) <https://doi.org/10.1016/j.ress.2021.107815>. Includes helpers
    to optimise toy models implemented in C++, visualise indices with
    uncertainty quantification, and derive reliability-oriented sensitivity
    measures based on failure probabilities.
    It is further detailed in Logosha, Maumy and Bertrand (2022) 
    <https://doi.org/10.1063/5.0246026> and (2023) <https://doi.org/10.1063/5.0246024> or in Bertrand, 
    Logosha and Maumy (2024) <https://hal.science/hal-05371803>, 
    <https://hal.science/hal-05371795> and <https://hal.science/hal-05371798>.
    
This site was created by F. Bertrand and the examples reproduced on it were created by F. Bertrand, E. Logosha and M. Maumy.

## Installation

You can install the latest version of the Sobol4R package from [github](https://github.com) with:


``` r
devtools::install_github("fbertran/Sobol4R")
```

## Two complementary analysis paths

Sobol4R exposes two ways to compute sensitivity indices depending on your
workflow:

* **Reuse the `sensitivity` package estimators.** Provide pre-built designs and
  call `sensitivity::sobol()` or `sensitivity::sobol2007()`; the `autoplot()`
  methods in Sobol4R will visualise those results without changing your
  existing code.
* **Use the in-package Saltelli re-implementation.** The `sobol_design()` and
  `sobol_indices()` helpers build the matrices, run the model, and return
  a `sobol_result` object that can be summarised or plotted directly, with
  optional bootstrap quantiles for noisy simulations.

The README examples below demonstrate the second path, while the earlier
"Context and non random case" section illustrates interoperability with
`sensitivity`.

## Motivation

The methodology implemented in **Sobol4R** builds upon the stochastic Sobol
analysis described by Lebrun et al. (2021) in *Reliability Engineering &
System Safety*. The paper proposes to combine replicated simulator
runs with Sobol estimators to account for intrinsic noise. The package mirrors
this workflow:

1. Generate two Monte Carlo designs (A and B matrices).
2. Evaluate the stochastic simulator several times to estimate the mean
   response and the noise variance.
3. Derive first-order and total-order Sobol indices and visualise the results.

The package is also friendly with the `sensitivity` package and shows how to use the Sobol' indices estimators provided in this package to increase the capabilities of this `Sobol4R` package.

## Basic usage


``` r
library(Sobol4R)
set.seed(123)
design <- sobol_design(n = 256, d = 3, lower = rep(-pi, 3), upper = rep(pi, 3),
                       quasi = TRUE)
result <- sobol_indices(ishigami_model, design, replicates = 4,
                        keep_samples = TRUE)
result$data
```

The resulting object stores the Monte Carlo variance estimate, the average
noise variance across replications, and the Sobol indices. Diagnostic plots are
available through the provided `autoplot()` method:


``` r
autoplot(result)
```

When `keep_samples = TRUE`, bootstrap resamples quantify the estimator
uncertainty. The helper `summarise_sobol()` produces tidy quantiles that can be
visualised directly:


``` r
autoplot(result, show_uncertainty = TRUE, probs = c(0.1, 0.9), bootstrap = 100)
```

## Reliability metrics

The paper highlights the need to quantify failure probabilities associated with
critical performance levels. The package exposes a helper for that task:


``` r
set.seed(321)
simulated <- ishigami_model(matrix(runif(3000, -pi, pi), ncol = 3))
estimate_failure_probability(simulated, threshold = -1)
```

When the simulator is stochastic and `sobol_indices()` stores the replicated
samples (`keep_samples = TRUE`), the same Monte Carlo budget can be recycled to
derive failure-indicator Sobol indices:


``` r
failure <- sobol_reliability(result, threshold = -1)
failure$failure_probability
autoplot(failure, show_uncertainty = TRUE, probs = c(0.1, 0.9), bootstrap = 200)
```

## Combined usage with the `sensitivity` package

### Test case : the non-monotonic Sobol g-function

The method of Sobol requires two samples. In this reference case there are eight variables, all following the uniform distribution on $[0,1]$.


``` r
if(require(sensitivity)){
n <- 50000
p <- 8
X1_1 <- data.frame(matrix(runif(p * n), nrow = n))
X2_1 <- data.frame(matrix(runif(p * n), nrow = n))
}
```


``` r
if(require(sensitivity)){
set.seed(4669)
gensol1 <- sobol4r_design(
  X1    = X1_1,
  X2    = X2_1,
  order = 2,
  nboot = 100
)

Y1 <- sobol_g_function(gensol1$X)
x1 <- sensitivity::tell(gensol1, Y1)
print(x1)
}
```


``` r
if(require(sensitivity)){
autoplot(x1, ncol = 1)
}
```



You can also use the `sobol_example_g_deterministic()` wrapper for this example.


``` r
if(require(sensitivity)){ex1_results <- sobol_example_g_deterministic()
print(ex1_results)
}
```


``` r
if(require(sensitivity)){
autoplot(ex1_results, ncol = 1)
}
```



### Vignettes

There are more insights and examples in the vignettes.


``` r
vignette("Sobol_RV_five_examples", package = "Sobol4R")
vignette("Sobol4R_vignette_stochastic", package = "Sobol4R")
vignette("Sobol4R_vignette_process", package = "Sobol4R")
vignette("simmer_MM1_Sobol_example", package = "Sobol4R")
```



