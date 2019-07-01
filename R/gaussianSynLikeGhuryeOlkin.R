#' Estimating the Gaussian synthetic likelihood with an unbiased estimator
#'
#' @description This function computes an unbiased, nonnegative estimate of a normal density 
#' function from simulations assumed to be drawn from it. See Price et al. (2018) and Ghurye 
#' and Olkin (1969).
#'
#' @inheritParams gaussianSynLike
#'
#' @return             The estimated synthetic (log) likelihood value.
#'
#' @references
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2017.1302882}
#'
#' Ghurye, S. G., & Olkin, I. (1969).
#' Unbiased estimation of some multivariate probability densities and related functions.
#' The Annals of Mathematical Statistics. \url{https://projecteuclid.org/euclid.aoms/1177697501}
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
#' # unbiased estimate of the Gaussian synthetic likelihood
#' # (the likelihood estimator used in uBSL)
#' gaussianSynLikeGhuryeOlkin(y, x)
#' 
#' @seealso    \code{\link{gaussianSynLike}} for the standard synthetic likelihood estimator, 
#' \code{\link{semiparaKernelEstimate}} for the semi-parametric likelihood estimator.
#' @export
gaussianSynLikeGhuryeOlkin <- function(ssy, ssx, log = TRUE, verbose = FALSE) {
    d <- length(ssy)
	n <- nrow(ssx)
    mu <- colMeans(ssx)
    Sigma <- cov(ssx)
	psi <- (n-1) * Sigma - (ssy-mu) %*% t(ssy-mu) / (1-1/n)

	temp <- try(chol(psi))
	if (inherits(temp, 'try-error')) {
	    if (verbose) {
            cat('*** reject (cov(ssx) is not positive definite) ***\n')
        }
		loglike <- -Inf
	} else {
	    A <- wcon(d, n-2) - wcon(d, n-1) - 0.5*d*log(1-1/n)
	    B <- -0.5 * (n-d-2) * (log(n-1) + logdet(Sigma))
	    C <- 0.5 * (n-d-3) * logdet(psi)
	    loglike <- -0.5*d*log(2*pi) + A + B + C
	}
	if (!log) {
	    loglike <- exp(loglike)
	}

    return (loglike)
}

# log of c(k,nu) from Ghurye & Olkin (1969)
wcon <- function(k, nu) {
    cc <- -k*nu/2*log(2) - k*(k-1)/4*log(pi) - sum(lgamma(0.5*(nu-(1:k)+1)))
	cc
}

# calculating the log of the determinant
logdet <- function(A) {
    L <- chol(A)
	y <- 2 * sum(log(diag(L)))
	y
}
