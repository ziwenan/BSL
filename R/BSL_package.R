#' Bayesian synthetic likelihood
#' 
#' @description 
#' Bayesian synthetic likelihood (BSL, Price et al (2018)) is an alternative to standard,
#' non-parametric approximate Bayesian computation (ABC). BSL assumes a multivariate normal distribution
#' for the summary statistic likelihood and it is suitable when the distribution of the model summary
#' statistics is sufficiently regular.
#' 
#' In this package, a Metropolis Hastings Markov chain Monte Carlo (MH-MCMC) implementation of BSL is available.
#' We also include an implementation of BSLasso (An et al., 2018a) which is computationally more efficient when the
#' dimension of the summary statistic is high, and semiBSL (An et al. 2018b) which relaxes the normality assumption
#' to an extent and maintain the computational advantages of BSL without any tuning. 
#'
#' Parallel computing is supported through the \code{foreach} package and users can specify their own parallel
#' backend by using packages like \code{doParallel} or \code{doMC}.
#'
#' The main functionality is available through
#'
#' \itemize{
#'      \item \code{\link{bsl}}: The general function to perform BSL, BSLasso and semiBSL (with or without parallel computing).
#'      \item \code{\link{selectPenalty}}: A function to select the penalty for BSLasso.
#'}
#'
#' Several examples have also been included. These examples can be used to reproduce the results of An et al. (2018a).
#'
#' \itemize{
#'      \item \code{\link{ma2}}: The MA(2) example from An et al. (2018a).
#'      \item \code{\link{mgnk}}: The multivariate G&K example from An et al. (2018a).
#'      \item \code{\link{cell}}: The cell biology example from Price et al. (2018) and An et al. (2018a).
#'}
#'
#' Extensions to this package are planned. 
#'
#' @references
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2017.1302882}
#'
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018a). Accelerating Bayesian synthetic 
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#'
#' An, Z., Nott, D. J. &  Drovandi, C. (2018b). Robust Bayesian Synthetic Likelihood via
#' a Semi-Parametric Approach. ArXiv Preprint \url{https://arxiv.org/abs/1809.05800}
#'
#' @author Ziwen An, Leah F. South and Christopher C. Drovandi
"_PACKAGE"
#> [1] "_PACKAGE"