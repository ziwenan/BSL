#' Estimating the Gaussian synthetic likelihood
#'
#' @description This function estimates the Gaussian synthetic likelihood function of Wood (2010).
#' Shrinkage on the Gaussian covariance matrix is also available (see An et al 2018).
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
#' @return             The estimated (log) likelihood value.
#'
#' @references
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2017.1302882}
#'
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating Bayesian synthetic 
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#'
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://doi.org/10.1198/016214508000000021}
#'
#' @examples
#' data(ma2)
#' y <- ma2$data # the observed data
#'
#' # function that simulates an ma2 time series
#' simulate_ma2 <- function(theta, L = 50) {
#'     rand <- rnorm(L + 2)
#'     y <- rand[3 : (L+2)] + theta[1] * rand[2 : (L+1)] + theta[2] * rand[1 : L]
#'     return(y)
#' }
#'
#' theta_true <- c(0.6, 0.2)
#' x <- matrix(0, 300, 50)
#' set.seed(100)
#' for(i in 1:300) x[i, ] <- simulate_ma2(theta_true)
#'
#' # the standard Gaussian synthetic likelihood (the likelihood estimator used in BSL)
#' gaussianSynLike(y, x)
#' # the Gaussian synthetic likelihood with glasso shrinkage estimation (the likelihood estimator used in BSLasso)
#' gaussianSynLike(y, x, shrinkage = 'glasso', penalty = 0.1)
#'
#' @export
gaussianSynLike <- function(ssy, ssx, shrinkage = NULL, penalty = NULL, standardise = FALSE, log = TRUE, verbose = FALSE) {
    if (is.null(shrinkage) && !is.null(penalty)) {
        warning('"penalty" will not ignored since no shrinkage method is specified')
    }
    if (!is.null(shrinkage) && is.null(penalty)) {
        stop('"penalty" must be specified to provoke shrinkage method')
    }
    if (is.null(penalty) & standardise) {
        warning('standardisation is only supported when shrinkage is "glasso"')
    }
    if (is.null(shrinkage)) { # BSL if no shrinkage
        mu <- colMeans(ssx)
        Sigma <- cov(ssx)
    } else { # BSL with shrinkage (glasso or Warton)
        mu <- colMeans(ssx)
        if (shrinkage == 'glasso') {
            if (!standardise) { # use graphical lasso without standardisation
                S <- cov(ssx)
                gl <- glasso(S, rho = penalty)
                Sigma <- gl$w
            } else { # standardise the summary statistics before passing into the graphical lasso function
                n <- nrow(ssx)
                ns <- ncol(ssx)
                std <- apply(ssx, MARGIN = 2, FUN = sd)
                ssx_std <- (ssx - matrix(mu, n, ns, byrow = TRUE)) / matrix(std, n, ns, byrow = TRUE)
                S <- cov(ssx_std)
                gl <- glasso(S, rho = penalty, penalize.diagonal = FALSE) # do not penalise the diagonal entries since we want the correlation matrix
                corr <- gl$w
                Sigma <- outer(std, std) * corr
            }
        } else if (shrinkage == 'Warton') {
            S <- cov(ssx)
            Sigma <- covWarton(S, penalty)
        } else {
            stop('shrinkage must be one of "glasso" and "Warton"')
        }
    }
    loglike <- try(mvtnorm::dmvnorm(ssy, mean = mu, sigma = Sigma, log = log))
    if (inherits(loglike, 'try-error')) {
        if (verbose) {
            cat('*** reject (probably singular cov(ssx) matrix) ***\n')
        }
        return (-Inf)
    }
    return (loglike)
}
