#' Estimating the semi-parametric joint likelihood
#'
#' @description This function computes the semi-parametric likelihood estimator of An et al (2018).
#' Kernel density estimates are used for modelling each univariate marginal distribution, and the
#' dependence structure between summaries are captured using a Gaussian copula.
#'
#' @param ssy          The observed summary statistic.
#' @param ssx          A matrix of the simulated summary statistics. The number of rows is the same
#' as the number of simulations per iteration.
#' @param kernel       A string argument indicating the smoothing kernel to pass into
#' \code{density} for estimating the marginal distribution of each summary statistic. Only ``gaussian"
#' and ``epanechnikov" are available. The default is ``gaussian".
#' @param shrinkage     A string argument indicating which shrinkage method to be used on the
#' correlation matrix of the Gaussian copula. The default is \code{NULL}, which means no shrinkage
#' is used. Current options are ``glasso" for graphical lasso and ``Warton" for the ridge
#' regularisation method of Warton (2008).
#' @inheritParams bsl
#' @param log          A logical argument indicating if the log of the likelihood is given as the result.
#' The default is \code{TRUE}.
#'
#' @return             The estimated synthetic (log) likelihood value.
#'
#' @references
#' An, Z., Nott, D. J. &  Drovandi, C. (2018). Robust Bayesian Synthetic Likelihood via
#' a Semi-Parametric Approach. \url{https://arxiv.org/abs/1809.05800}
#'
#' Friedman, J., Hastie, T., Tibshirani, R. (2008). Sparse inverse covariance estimation with
#' the graphical lasso. Biostatistics. \url{https://doi.org/10.1093/biostatistics/kxm045}
#'
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://www.tandfonline.com/doi/abs/10.1198/016214508000000021}
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
#' # the default semi-parametric synthetic likelihood estimator of semiBSL
#' semiparaKernelEstimate(y, x)
#' # using shrinkage on the correlation matrix of the Gaussian copula is also possible
#' semiparaKernelEstimate(y, x, shrinkage = 'Warton', penalty = 0.6)
#'
#' @seealso    \code{\link{gaussianSynLike}} for the standard synthetic likelihood estimator, 
#' \code{\link{gaussianSynLikeGhuryeOlkin}} for the unbiased synthetic likelihood estimator.
#' @export
semiparaKernelEstimate <- function (ssy, ssx, kernel = 'gaussian', shrinkage = NULL, penalty = NULL, log = TRUE) {
    if (is.null(shrinkage) && !is.null(penalty)) {
        warning('"penalty" will not ignored since no shrinkage method is specified')
    }
    if (!is.null(shrinkage) && is.null(penalty)) {
        stop('"penalty" must be specified to provoke shrinkage method')
    }

    n <- nrow(ssx)
    ns <- ncol(ssx)

    pdfy <- yu <- numeric(ns)
    for (j in 1 : ns) {
        f <- density(ssx[, j], kernel = kernel, n = 512, from = min(ssx[,j],ssy[j]), to = max(ssx[,j],ssy[j]))
        approxy <- approx(f$x, f$y, ssy[j])
        pdfy[j] <- approxy$y
        yu[j] <- mean(kernelCDF((ssy[j] - ssx[, j]) / f$bw, kernel))
    }

    # Gaussian rank correlation
	rhohat <- gaussianRankCorr(ssx, TRUE)

    if (!is.null(shrinkage)) {
        RHOHAT <- p2P(rhohat)
        Sigma <- switch(shrinkage,
                        'glasso' = glasso(RHOHAT, rho = penalty, penalize.diagonal = FALSE)$w,
                        'Warton' = corrWarton(RHOHAT, penalty))
        rhohat <- P2p(Sigma)
    }

    # density
    copula <- normalCopula(rhohat, dim = ns, dispstr = 'un')
	if (log) {
	    f <- dCopula(yu, copula, log = TRUE) + sum(log(pdfy))
	} else {
	    f <- dCopula(yu, copula, log = FALSE) * prod(pdfy)
	}

    return(f)
}
