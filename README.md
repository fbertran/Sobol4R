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
#>   parameter first_order total_order
#> 1        X1           0           0
#> 2        X2           0           0
#> 3        X3           0           0
```

The resulting object stores the Monte Carlo variance estimate, the average
noise variance across replications, and the Sobol indices. Diagnostic plots are
available through the provided `autoplot()` method:


``` r
autoplot(result)
#> Error in autoplot(result): could not find function "autoplot"
```

When `keep_samples = TRUE`, bootstrap resamples quantify the estimator
uncertainty. The helper `summarise_sobol()` produces tidy quantiles that can be
visualised directly:


``` r
autoplot(result, show_uncertainty = TRUE, probs = c(0.1, 0.9), bootstrap = 100)
#> Error in autoplot(result, show_uncertainty = TRUE, probs = c(0.1, 0.9), : could not find function "autoplot"
```

## Reliability metrics

The paper highlights the need to quantify failure probabilities associated with
critical performance levels. The package exposes a helper for that task:


``` r
set.seed(321)
simulated <- ishigami_model(matrix(runif(3000, -pi, pi), ncol = 3))
estimate_failure_probability(simulated, threshold = -1)
#> $probability
#> [1] 0.087
#> 
#> $variance
#> [1] 7.9431e-05
```

When the simulator is stochastic and `sobol_indices()` stores the replicated
samples (`keep_samples = TRUE`), the same Monte Carlo budget can be recycled to
derive failure-indicator Sobol indices:


``` r
failure <- sobol_reliability(result, threshold = -1)
failure$failure_probability
#> [1] 0.08984375
autoplot(failure, show_uncertainty = TRUE, probs = c(0.1, 0.9), bootstrap = 200)
#> Error in autoplot(failure, show_uncertainty = TRUE, probs = c(0.1, 0.9), : could not find function "autoplot"
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
#> 
#> Call:
#> sensitivity::sobol(model = NULL, X1 = X1, X2 = X2, order = order,     nboot = nboot)
#> 
#> Model runs: 1850000 
#> 
#> Sobol indices
#>            original          bias  std. error
#> X1     0.7196152898  0.0013178498 0.007451246
#> X2     0.1853898383 -0.0001958423 0.009955315
#> X3     0.0315602561  0.0003865755 0.010215090
#> X4     0.0165983351  0.0006435924 0.009825790
#> X5     0.0071457842  0.0005218573 0.009914644
#> X6     0.0070476581  0.0005221334 0.009899108
#> X7     0.0072387630  0.0005036510 0.009885607
#> X8     0.0069760096  0.0004883120 0.009909057
#> X1*X2  0.0482402280 -0.0008197148 0.011798374
#> X1*X3  0.0003008437 -0.0006354620 0.010346606
#> X1*X4 -0.0054547860 -0.0005409433 0.010174239
#> X1*X5 -0.0069436347 -0.0005097920 0.009898335
#> X1*X6 -0.0070049660 -0.0005216520 0.009901355
#> X1*X7 -0.0069548448 -0.0005165511 0.009891716
#> X1*X8 -0.0068674169 -0.0005250381 0.009898584
#> X2*X3 -0.0058530637 -0.0006812381 0.009961972
#> X2*X4 -0.0059294545 -0.0004713903 0.009984341
#> X2*X5 -0.0070401939 -0.0005053830 0.009910392
#> X2*X6 -0.0070285618 -0.0005149528 0.009899016
#> X2*X7 -0.0069326967 -0.0005195100 0.009906751
#> X2*X8 -0.0070601319 -0.0005246410 0.009903606
#> X3*X4 -0.0067304606 -0.0005017321 0.009898065
#> X3*X5 -0.0070185807 -0.0005154295 0.009901753
#> X3*X6 -0.0070280907 -0.0005135826 0.009903321
#> X3*X7 -0.0070005474 -0.0005151754 0.009900567
#> X3*X8 -0.0069922906 -0.0005170020 0.009901918
#> X4*X5 -0.0070143438 -0.0005165487 0.009902715
#> X4*X6 -0.0070074685 -0.0005170955 0.009902002
#> X4*X7 -0.0070098377 -0.0005165402 0.009902303
#> X4*X8 -0.0070059467 -0.0005179379 0.009903008
#> X5*X6 -0.0070067930 -0.0005160748 0.009902998
#> X5*X7 -0.0070059994 -0.0005161184 0.009903020
#> X5*X8 -0.0070063914 -0.0005160893 0.009902842
#> X6*X7 -0.0070065347 -0.0005160524 0.009902864
#> X6*X8 -0.0070060415 -0.0005159716 0.009902961
#> X7*X8 -0.0070048263 -0.0005159800 0.009902875
#>          min. c.i.  max. c.i.
#> X1     0.700120385 0.73334360
#> X2     0.163759252 0.20599792
#> X3     0.008547491 0.04573165
#> X4    -0.006227570 0.03128409
#> X5    -0.015044519 0.02224057
#> X6    -0.015191431 0.02197302
#> X7    -0.014887579 0.02219752
#> X8    -0.015196240 0.02193597
#> X1*X2  0.026200069 0.07361523
#> X1*X3 -0.014125796 0.02290868
#> X1*X4 -0.021045640 0.01734604
#> X1*X5 -0.021899769 0.01530898
#> X1*X6 -0.022033039 0.01534886
#> X1*X7 -0.021889386 0.01535810
#> X1*X8 -0.021701484 0.01536219
#> X2*X3 -0.021921145 0.01612148
#> X2*X4 -0.021220929 0.01626216
#> X2*X5 -0.022020473 0.01527561
#> X2*X6 -0.021921481 0.01526160
#> X2*X7 -0.021840665 0.01539284
#> X2*X8 -0.021959503 0.01527697
#> X3*X4 -0.021619338 0.01561685
#> X3*X5 -0.021932111 0.01530112
#> X3*X6 -0.021965128 0.01526305
#> X3*X7 -0.021926020 0.01528115
#> X3*X8 -0.021890026 0.01530904
#> X4*X5 -0.021939624 0.01529773
#> X4*X6 -0.021924396 0.01530570
#> X4*X7 -0.021937676 0.01529542
#> X4*X8 -0.021931999 0.01529803
#> X5*X6 -0.021927352 0.01530282
#> X5*X7 -0.021926839 0.01530200
#> X5*X8 -0.021928525 0.01530312
#> X6*X7 -0.021928273 0.01530131
#> X6*X8 -0.021927491 0.01530197
#> X7*X8 -0.021926479 0.01530484
```


``` r
if(require(sensitivity)){
autoplot(x1, ncol = 1)
}
#> Error in autoplot(x1, ncol = 1): could not find function "autoplot"
```



You can also use the `sobol_example_g_deterministic()` wrapper for this example.


``` r
if(require(sensitivity)){ex1_results <- sobol_example_g_deterministic()
print(ex1_results)
}
#> 
#> Call:
#> sensitivity::sobol(model = NULL, X1 = X1, X2 = X2, order = order,     nboot = nboot)
#> 
#> Model runs: 1850000 
#> 
#> Sobol indices
#>            original          bias  std. error
#> X1     0.7245997507  1.318649e-04 0.006865099
#> X2     0.1852412158 -6.379462e-04 0.009725422
#> X3     0.0321041221 -3.943572e-04 0.009939738
#> X4     0.0150373622 -3.716233e-04 0.009571601
#> X5     0.0073639355 -5.240577e-04 0.009690646
#> X6     0.0073304377 -5.140176e-04 0.009697496
#> X7     0.0072934310 -5.369366e-04 0.009679297
#> X8     0.0070625492 -5.292390e-04 0.009661789
#> X1*X2  0.0459216617  8.939108e-05 0.010932749
#> X1*X3 -0.0006600465  6.814819e-04 0.010010933
#> X1*X4 -0.0056037444  4.901684e-04 0.009860488
#> X1*X5 -0.0070363484  5.187023e-04 0.009676301
#> X1*X6 -0.0071411552  5.319393e-04 0.009690812
#> X1*X7 -0.0072518362  5.303046e-04 0.009672163
#> X1*X8 -0.0070777721  5.186929e-04 0.009676468
#> X2*X3 -0.0051274794  5.279125e-04 0.009702252
#> X2*X4 -0.0060860874  5.210190e-04 0.009681757
#> X2*X5 -0.0071063957  5.147476e-04 0.009680715
#> X2*X6 -0.0071163219  5.256144e-04 0.009679855
#> X2*X7 -0.0070620281  5.299198e-04 0.009678650
#> X2*X8 -0.0071567767  5.166541e-04 0.009686683
#> X3*X4 -0.0073116996  5.314332e-04 0.009671714
#> X3*X5 -0.0071206761  5.265249e-04 0.009682280
#> X3*X6 -0.0071350887  5.241734e-04 0.009680955
#> X3*X7 -0.0071632203  5.264482e-04 0.009681046
#> X3*X8 -0.0071109279  5.241054e-04 0.009682418
#> X4*X5 -0.0071437469  5.248274e-04 0.009679986
#> X4*X6 -0.0071379129  5.276560e-04 0.009680886
#> X4*X7 -0.0071596546  5.255998e-04 0.009681341
#> X4*X8 -0.0071300368  5.260610e-04 0.009682262
#> X5*X6 -0.0071348129  5.263360e-04 0.009681340
#> X5*X7 -0.0071382804  5.262539e-04 0.009681217
#> X5*X8 -0.0071340327  5.262561e-04 0.009681110
#> X6*X7 -0.0071357204  5.261516e-04 0.009681295
#> X6*X8 -0.0071339651  5.264348e-04 0.009681123
#> X7*X8 -0.0071370385  5.263348e-04 0.009681299
#>          min. c.i.  max. c.i.
#> X1     0.711583661 0.73855259
#> X2     0.163919891 0.20792418
#> X3     0.012874359 0.05265936
#> X4    -0.002765501 0.03471590
#> X5    -0.010879190 0.02837068
#> X6    -0.010964117 0.02838793
#> X7    -0.010989042 0.02830642
#> X8    -0.010934969 0.02789333
#> X1*X2  0.026094148 0.06869013
#> X1*X3 -0.022980303 0.01844932
#> X1*X4 -0.026644889 0.01263961
#> X1*X5 -0.027999214 0.01120229
#> X1*X6 -0.028192638 0.01110205
#> X1*X7 -0.028265971 0.01093724
#> X1*X8 -0.028051234 0.01119551
#> X2*X3 -0.026495177 0.01365939
#> X2*X4 -0.027207071 0.01195456
#> X2*X5 -0.028066083 0.01115791
#> X2*X6 -0.028098511 0.01109647
#> X2*X7 -0.028017158 0.01116083
#> X2*X8 -0.028157933 0.01108562
#> X3*X4 -0.028275534 0.01102687
#> X3*X5 -0.028100457 0.01114042
#> X3*X6 -0.028110685 0.01112504
#> X3*X7 -0.028151981 0.01111060
#> X3*X8 -0.028109093 0.01117038
#> X4*X5 -0.028127750 0.01111732
#> X4*X6 -0.028127090 0.01112126
#> X4*X7 -0.028150168 0.01109302
#> X4*X8 -0.028120750 0.01114550
#> X5*X6 -0.028121862 0.01112445
#> X5*X7 -0.028124558 0.01111832
#> X5*X8 -0.028120227 0.01112288
#> X6*X7 -0.028121656 0.01112018
#> X6*X8 -0.028121004 0.01112360
#> X7*X8 -0.028124733 0.01112035
```


``` r
if(require(sensitivity)){
autoplot(ex1_results, ncol = 1)
}
#> Error in autoplot(ex1_results, ncol = 1): could not find function "autoplot"
```



### Vignettes

There are more insights and examples in the vignettes.


``` r
vignette("Sobol_RV_five_examples", package = "Sobol4R")
vignette("Sobol4R_vignette_stochastic", package = "Sobol4R")
vignette("Sobol4R_vignette_process", package = "Sobol4R")
vignette("simmer_MM1_Sobol_example", package = "Sobol4R")
```



