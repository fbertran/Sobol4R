#' Sobol4R-package
#'
#'
#'
#' @references  Elizaveta Logosha, Myriam Maumy, Frederic Bertrand; Confidence interval determination using discrete event simulations for real estate sales case. AIP Conf. Proc. 31 March 2025; 3182 (1): 100008. \doi{10.1063/5.0246026}.
#' 
#' Elizaveta Logosha, Myriam Maumy, Frédéric Bertrand; Sensitivity analysis of stochastic simulator in the case of sales date prediction. AIP Conf. Proc. 31 March 2025; 3182 (1): 100001. \doi{10.1063/5.0246024}.
#' 
#' Frédéric Bertrand, Elizaveta Logosha, Myriam Maumy-Bertrand. Extension of sensitivity analysis to uncertainties in distribution parameters. 32nd Conference on Intelligent Systems for Molecular Biology, International Society for Computational Biology, Jul 2024, Montreal (QC), Canada. \url{https://hal.science/hal-05371795}.
#' 
#' Frédéric Bertrand, Elizaveta Logosha, Myriam Maumy-Bertrand. Sobol4RV: Global Sensitivity Analysis in Several Random Settings. BioC 2024, BioConductor, Jul 2024, Grand Rapids, MI, United States. \url{https://hal.science/hal-05371803}
#' 
#' Frédéric Bertrand, Elizaveta Logosha, Myriam Maumy-Bertrand. Global Sensitivity Analysis in Several Random Settings. 2024 Joint Statistical Meetings, American Statistical Association, Aug 2024, Portland (OR), United States. \url{https://hal.science/hal-05371798}.
#'
#' @seealso [sobol4r_design()], [sobol4r_qoi_indices()], 
#'   `vignette("Sobol_RV_five_examples", package = "Sobol4R")`,
#'   `vignette("Sobol4R_vignette_stochastic", package = "Sobol4R")`,
#'   `vignette("Sobol4R_vignette_process", package = "Sobol4R")` and
#'   `vignette("simmer_MM1_Sobol_example", package = "Sobol4R")`.
#' 
#' @examples
#' ex1_results <- sobol_example_g_deterministic(n=100, nboot=10) 
#' print(ex1_results)
#' autoplot(ex1_results, ncol = 1)
#' rm(ex1_results)
#' 
"_PACKAGE"

#' @importFrom rlang .data is_installed
#' @importFrom ggplot2 aes autoplot facet_wrap geom_col geom_errorbar geom_hline ggplot labs position_dodge theme theme_minimal
#' @importFrom graphics barplot
#' @useDynLib Sobol4R, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
