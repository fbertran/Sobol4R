#' Sobol4R-package
#'
#' A minimal benchmark for users with GPU access is:
#'
#' @examples
#' 
#' library(Sobol4R)
#' if (!requireNamespace("sensitivity", quietly = TRUE)) {
#' 
#' set.seed(123)
#' n <- 50000
#' X1 <- data.frame(matrix(runif(2 * n), nrow = n))
#' X2 <- data.frame(matrix(runif(2 * n), nrow = n))
#' 
#' res_det <- sobol4r_run(
#'   model = sobol4r_g2,
#'   X1 = X1,
#'   X2 = X2,
#'   order = 2,
#'   nboot = 100
#' )
#' 
#' print(res_det)
#' ggplot2::ggplot(res_det)
#' }
#'
#' @aliases Sobol4R-package Sobol4R NULL
#'
#' @references  Elizaveta Logosha, Myriam Maumy, Frederic Bertrand; Confidence interval determination using discrete event simulations for real estate sales case. AIP Conf. Proc. 31 March 2025; 3182 (1): 100008. <doi.org:10.1063/5.0246026>.
#' 
#' Elizaveta Logosha, Myriam Maumy, Frédéric Bertrand; Sensitivity analysis of stochastic simulator in the case of sales date prediction. AIP Conf. Proc. 31 March 2025; 3182 (1): 100001. <doi.org:10.1063/5.0246024>
#' 
#' Frédéric Bertrand, Elizaveta Logosha, Myriam Maumy-Bertrand. Extension of sensitivity analysis to uncertainties in distribution parameters. 32nd Conference on Intelligent Systems for Molecular Biology, International Society for Computational Biology, Jul 2024, Montreal (QC), Canada. <https://hal.science/hal-05371795>
#' 
#' Frédéric Bertrand, Elizaveta Logosha, Myriam Maumy-Bertrand. Sobol4RV: Global Sensitivity Analysis in Several Random Settings. BioC 2024, BioConductor, Jul 2024, Grand Rapids, MI, United States. <https://hal.science/hal-05371803>
#' 
#' Frédéric Bertrand, Elizaveta Logosha, Myriam Maumy-Bertrand. Global Sensitivity Analysis in Several Random Settings. 2024 Joint Statistical Meetings, American Statistical Association, Aug 2024, Portland (OR), United States. <https://hal.science/hal-05371798>
#'
#' @seealso TODO
#' 
"_PACKAGE"

#' @importFrom stats rbinom rexp rnorm runif
#' @importFrom sensitivity sobol tell
# #' @useDynLib Sobol4R, .registration = TRUE
# #' @importFrom Rcpp evalCpp
NULL


