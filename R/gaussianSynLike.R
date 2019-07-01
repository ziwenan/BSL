#' Estimating the Gaussian synthetic likelihood
#'
#' @description This function estimates the Gaussian synthetic likelihood function of Wood (2010).
#' Shrinkage on the Gaussian covariance matrix is also available (see An et al 2019).
#'
#' @param ssy          The observed summary statisic.
#' @param ssx          A matrix of the simulated summary statistics. The number of rows is the same
#' as the number of simulations per iteration.
#' @param log          A logical argument indicating if the log of likelihood is given as the result.
#' The default is \code{TRUE}.
#' @param verbose      A logical argument indicating whether an error message should be printed if
#' the function fails to compute a likelihood. The default is \code{FALSE}.
#' @inheritParams bsl
#'
#' @return             The estimated synthetic (log) likelihood value.
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
#' Friedman, J., Hastie, T., Tibshirani, R. (2008). Sparse inverse covariance estimation with
#' the graphical lasso. Biostatistics. \url{https://doi.org/10.1093/biostatistics/kxm045}
#'
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://doi.org/10.1198/016214508000000021}
#'
#' @examples
#' data(ma2)
#' y <- ma2$data # the observed data
#' 
#' theta_true <- c(0.6, 0.2)
#' x <- matrix(0, 300, 50)
#' set.seed(100)
#' for(i in 1:300) x[i, ] <- ma2_sim(theta_true, 50)
#' 
#' # the standard Gaussian synthetic likelihood (the likelihood estimator used in BSL)
#' gaussianSynLike(y, x)
#' # the Gaussian synthetic likelihood with glasso shrinkage estimation
#' # (the likelihood estimator used in BSLasso)
#' gaussianSynLike(y, x, shrinkage = 'glasso', penalty = 0.1)
#' # the Gaussian synthetic likelihood with Warton shrinkage estimation
#' gaussianSynLike(y, x, shrinkage = 'Warton', penalty = 0.9)
#'
#' @seealso    \code{\link{gaussianSynLikeGhuryeOlkin}} for the unbiased synthetic likelihood estimator, 
#' \code{\link{semiparaKernelEstimate}} for the semi-parametric likelihood estimator.
#' @export
gaussianSynLike <- function(ssy, ssx, shrinkage = NULL, penalty = NULL, standardise = FALSE, GRC = FALSE, log = TRUE, verbose = FALSE) {
    if (is.null(shrinkage) && !is.null(penalty)) {
        warning('"penalty" will be ignored since no shrinkage method is specified')
    }
    if (!is.null(shrinkage) && is.null(penalty)) {
        stop('"penalty" must be specified to provoke shrinkage method')
    }
    if (shrinkage != 'glasso' && standardise) {
        warning("standardisation is only supported in BSLasso")
    }
    if (is.null(shrinkage)) { # BSL if no shrinkage
        mu <- colMeans(ssx)
		if (GRC) {
		    std <- apply(ssx, MARGIN = 2, FUN = sd)
		    corr <- gaussianRankCorr(ssx)
			Sigma <- cor2cov(corr, std)
		} else {
		    Sigma <- cov(ssx)
		}
    } else { # BSL with shrinkage (glasso or Warton)
        mu <- colMeans(ssx)
        if (shrinkage == 'glasso') {
            if (!standardise) { # use graphical lasso without standardisation
			    if (GRC) {
				    std <- apply(ssx, MARGIN = 2, FUN = sd)
		            corr <- gaussianRankCorr(ssx)
			        S <- cor2cov(corr, std)
				} else {
				    S <- cov(ssx)
				}
                gl <- glasso(S, rho = penalty)
                Sigma <- gl$w
            } else { # standardise the summary statistics before passing into the graphical lasso function
                n <- nrow(ssx)
                ns <- ncol(ssx)
                std <- apply(ssx, MARGIN = 2, FUN = sd)
                ssx_std <- (ssx - matrix(mu, n, ns, byrow = TRUE)) / matrix(std, n, ns, byrow = TRUE)
				if (GRC) {
		            corr <- gaussianRankCorr(ssx_std)
			        S <- cor2cov(corr, apply(ssx_std, MARGIN = 2, FUN = sd))
				} else {
				    S <- cov(ssx_std)
				}
                gl <- glasso(S, rho = penalty, penalize.diagonal = FALSE) # do not penalise the diagonal entries since we want the correlation matrix
                corr <- gl$w
                Sigma <- outer(std, std) * corr
            }
        } else if (shrinkage == 'Warton') {
		    if (GRC) {
			    std <- apply(ssx, MARGIN = 2, FUN = sd)
		        corr <- gaussianRankCorr(ssx)
			    S <- cor2cov(corr, std)
			} else {
			    S <- cov(ssx)
			}
            Sigma <- covWarton(S, penalty)
        } else {
            stop('shrinkage must be one of "glasso" and "Warton"')
        }
    }
    loglike <- try(mvtnorm::dmvnorm(ssy, mean = mu, sigma = Sigma, log = log))
    if (inherits(loglike, 'try-error')) {
        if (verbose) {
            cat('*** reject (probably singular covariance matrix) ***\n')
        }
        return (-Inf)
    }
    return (loglike)
}
