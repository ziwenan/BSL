#' Bayesian synthetic likelihood
#' 
#' @description 
#' Bayesian synthetic likelihood (BSL, Price et al (2018)) is an alternative to standard,
#' non-parametric approximate Bayesian computation (ABC). BSL assumes a multivariate normal distribution
#' for the summary statistic likelihood and it is suitable when the distribution of the model summary
#' statistics is sufficiently regular.
#' 
#' In this package, a Metropolis Hastings Markov chain Monte Carlo (MH-MCMC) implementation of BSL is available.
#' We also include implementations of three methods (BSL, uBSL and semiBSL) and two shrinkage estimations 
#' (graphical lasso and Warton's estimation).
#'
#' Methods:
#' (1) BSL (Price et al., 2018), which is the standard form of Bayesian synthetic likelihood, assumes the summary 
#' statistic is roughly multivariate normal; (2) uBSL (Price et al., 2018), which uses an unbiased estimator to 
#' the normal density; and (3) semiBSL (An et al. 2018), which relaxes the normality assumption to an extent and 
#' maintains the computational advantages of BSL without any tuning. 
#'
#' Shrinkage estimations are designed particularly to reduce the number of simulations if method is BSL or semiBSL:
#' (1) graphical lasso (Friedman et al., 2008) finds a sparse precision matrix with a L1-regularised log-likelihood.
#' (An et al., 2019) use graphical lasso within BSL to bring down the number of simulations significantly when 
#' the dimension of the summary statistic is high; and (2) Warton's estimator (Warton 2008) penalises the 
#' correlation matrix and is straightforward to compute.
#'
#' Parallel computing is supported through the \code{foreach} package and users can specify their own parallel
#' backend by using packages like \code{doParallel} or \code{doMC}. Alternatively a vectorised simulation function
#' that simultaneously generates \code{n} simulation results is also supported.
#'
#' The main functionality is available through
#'
#' \itemize{
#'      \item \code{\link{bsl}}: The general function to perform BSL, uBSL, and semiBSL (with or without 
#'          parallel computing).
#'      \item \code{\link{selectPenalty}}: A function to select the penalty for BSL and semiBSL.
#'}
#'
#' Several examples have also been included. These examples can be used to reproduce the results of An et al. (2019), 
#' and can help practitioners learn how to use the package.
#'
#' \itemize{
#'      \item \code{\link{ma2}}: The MA(2) example from An et al. (2019).
#'      \item \code{\link{mgnk}}: The multivariate G&K example from An et al. (2019).
#'      \item \code{\link{cell}}: The cell biology example from Price et al. (2018) and An et al. (2019).
#'}
#'
#' Extensions to this package are planned. 
#'
#' @references
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2017.1302882}
#'
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2019). Accelerating Bayesian synthetic 
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#'
#' An, Z., Nott, D. J. &  Drovandi, C. (2018). Robust Bayesian Synthetic Likelihood via
#' a Semi-Parametric Approach. ArXiv Preprint \url{https://arxiv.org/abs/1809.05800}
#'
#' Friedman, J., Hastie, T., Tibshirani, R. (2008). Sparse inverse covariance estimation with
#' the graphical lasso. Biostatistics. \url{https://doi.org/10.1093/biostatistics/kxm045}
#'
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://doi.org/10.1198/016214508000000021}
#'
#' @author Ziwen An, Leah F. South and Christopher C. Drovandi
"_PACKAGE"
#> [1] "_PACKAGE"