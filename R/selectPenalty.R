#' Selecting BSLasso Penalty
#'
#' @description This is the main function for selecting the shrinkage (graphical lasso or Warton's estimation) penalty parameter
#' for method BSL or semiBSL based on a point estimate of the parameters. Parallel computing is supported with the R package 
#' \code{foreach}.
#'
#' @param ssy                   A summary statistic vector for the observed data.
#' @param n                     A vector of possible values of \code{n}, the number of simulations from the model per MCMC iteration for estimating the synthetic likelihood.
#' @param lambda_all            A list, with each entry containing the vector of penalty values to test for the corresponding choice of \code{n}.
#' @param theta                 A point estimate of the parameter value which all of the simulations will be based on.
#' @param M                     The number of repeats to use in estimating the standard deviation of the estimated log synthetic likelihood.
#' @param sigma                 The standard deviation of the log synthetic likelihood estimator to aim for, usually a value between 1 and 2. This reflects the mixing of a Markov chain.
#' @param method                A string argument indicating the method to be used. The default, ``BSL'', runs BSL.
#' ``semiBSL'' runs the semi-parametric BSL algorithm and is more robust to non-normal summary statistics.
#' @param shrinkage     A string argument indicating which shrinkage method to be used. Current options 
#' are ``glasso'' for the graphical lasso method of Friedman et al (2008) and ``Warton'' for the ridge regularisation method of Warton (2008).
#' @param parallelSim           A logical value indicating whether parallel computing should be used for simulation and summary statistic evaluation. Default is \code{FALSE}.
#' @param parallelSimArgs       A list of additional arguments to pass into the \code{foreach} function. Only used when parallelSim is \code{TRUE}, default is \code{NULL}.
#' @param parallelMain          A logical value indicating whether parallel computing should be used to computing the graphical lasso function. Default is \code{FALSE}.
#' @param verbose               A logical argument indicating whether the iteration numbers (\code{1:M}) should be printed to track progress. The default is \code{FALSE}.
#' @inheritParams bsl
#'
#' @return 				An object of class \code{penbsl} is returned, containing the following components:
#' \itemize{
#' \item \code{resultsDF}: A data frame containing the following:
#'    \itemize{
#'    \item \code{n}: The choices of \code{n} that were specified.
#'    \item \code{penalty}: The choices of the penalty that were specified.
#'    \item \code{sigma}: The standard deviation of the log synthetic likelihood estimator under the above choices.
#'    \item \code{sigmaOpt}: An indicator of whether it was the closest \code{sigma} to the desired one for each choice of \code{n}.
#'    }
#' \item \code{call}: The original code that was used to call the method.
#' }
#' The functions print() and plot() are both available for types of class \code{penbsl}.
#'
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2019). Accelerating Bayesian synthetic
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#' 
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://doi.org/10.1198/016214508000000021}
#'
#' @author    Ziwen An, Leah F. South and Christopher C. Drovandi
#' @seealso   \code{\link{ma2}}, \code{\link{cell}} and \code{\link{mgnk}} for examples. 
#' \code{\link{bsl}} for a function to run BSLasso after selecting the tuning parameter
#' and \code{\link{penbsl}} for functions related to visualisation.
#' @export
selectPenalty <- function(ssy, n, lambda_all, M, sigma, model, theta = model@theta0, 
    method = c('BSL', 'semiBSL')[1], shrinkage = c('glasso', 'Warton')[1], 
	standardise = FALSE, GRC = FALSE, parallelSim = FALSE, parallelSimArgs = NULL,  
	parallelMain = FALSE, verbose = TRUE, 
	fnSim, fnSum, simArgs, sumArgs) {
	if (standardise && (shrinkage != 'glasso' || method != 'BSL')) {
        warning('standardisation is only supported when method is "BSL" and shrinkage is "glasso"')
    }
	if (!parallelSim & !is.null(parallelSimArgs)) {
        warning("\"parallelSimArgs\" is omitted in serial computing")
    }
	if (!missing(fnSim) || !missing(fnSum) || !missing(simArgs) || !missing(sumArgs)) {
	    stop('fnSim, fnSum, simArgs and sumArgs are deprecated now, please use model instead')
	}
	stopifnot(class(model) == 'BSLMODEL')
	
    n <- as.vector(n)
    lambda_all <- as.list(lambda_all)
    N <- length(n)
    if (length(lambda_all) != N) {
        stop('lambda_all must be a list with the same length as n')
    }
    ns <- length(ssy)
    cl <- match.call()
	
	if (verbose) cat(paste0('*** selecting penalty with ', method, ' likelihood and ', shrinkage, ' shrinkage estimator ***\n'))
    n_max <- max(n)
    K <- max(sapply(lambda_all, length))
    logSL <- array(NA, c(M, N, K))

	# map the simulation function
	if (parallelSim) {
	    myFnSimSum <- function(n, theta) fn(model)$fnPar(n, theta, parallelSimArgs)
	} else {
	    myFnSimSum <- fn(model)$fn
	}

    for (m in 1 : M) {

        flush.console()
        if (verbose) cat('m =', m, '\n')

        # simulate with theta_prop and calculate summaries
		ssx <- myFnSimSum(n_max, theta)

        for (i in 1 : N) {
            n_curr <- n[i]
            ssx_curr <- ssx[sample(n_max, n_curr), ]

            if (!parallelMain) {
                for (k in 1 : K) {
				    lambda_curr = lambda_all[[i]][k]
				    if (is.na(lambda_curr)) {
				        next
				    }
					logSL[m, i, k] <- switch(method, 
					    'BSL' = gaussianSynLike(ssy, ssx_curr, shrinkage = shrinkage, penalty = lambda_curr, standardise = standardise, GRC = GRC, log = TRUE),
						'semiBSL' = semiparaKernelEstimate(ssy, ssx_curr, kernel = 'gaussian', shrinkage = shrinkage, penalty = lambda_curr, log = TRUE))
                }
            } else {
                logSL[m, i, ] <- foreach (k = 1 : K, .combine = rbind, .packages = 'glasso') %dopar% {
				    lambda_curr = lambda_all[[i]][k]
				    if (is.na(lambda_curr)) {
				        return(NA)
				    }
                    switch(method, 
					    'BSL' = gaussianSynLike(ssy, ssx_curr, shrinkage = shrinkage, penalty = lambda_curr, standardise = standardise, GRC = GRC, log = TRUE),
						'semiBSL' = semiparaKernelEstimate(ssy, ssx_curr, kernel = 'gaussian', shrinkage = shrinkage, penalty = lambda_curr, log = TRUE))
                }
            }

        }
    }

    resultsDF <- array(list(), N)
    for (i in 1 : N) {
        temp <- apply(logSL[, i, ], MARGIN = 2, FUN = sd)
        resultsDF[[i]] <- data.frame(n = n[i], penalty = lambda_all[[i]],
                              sigma = temp[!is.na(temp)], sigmaOpt = FALSE)
        resultsDF[[i]]$sigmaOpt[which.min(abs(resultsDF[[i]]$sigma - sigma))] <- TRUE
    }
    resultsDF <- do.call(rbind, resultsDF)

    results <- list(resultsDF = resultsDF, call = cl)
    class(results) <- 'penbsl'

    return(results)
}
