#' Performing BSL, BSLasso and semiBSL
#'
#' @description This is the main function for performing MCMC BSL, MCMC uBSL, MCMC BSLasso and MCMC semiBSL.
#' Parallel computing is supported with the R package \code{foreach}. Several input arguments in 2.0.0 or earlier versions
#' are replaced by \code{model} in the new update and will be deprecated in a future release.
#'
#' @param y				The observed data. Note this should be the raw dataset NOT the set of summary statistics.
#' @param n				The number of simulations from the model per MCMC iteration for estimating the synthetic likelihood.
#' @param M				The number of MCMC iterations.
#' @param model         A ``BSLMODEL'' object generated with function \code{BSLModel}. See \code{\link{BSLModel}}.
#' @param covRandWalk	The covariance matrix of a multivariate normal random walk proposal distribution used in the MCMC.
#' @param method        A string argument indicating the method to be used. The default, ``BSL'', runs standard BSL.
#' ``uBSL'' uses the unbiased estimator of a normal density of Ghurye and Olkin (1969). 
#' ``semiBSL'' runs the semi-parametric BSL algorithm and is more robust to non-normal summary statistics.
#' @param shrinkage     A string argument indicating which shrinkage method to be used. The default is \code{NULL},
#' which means no shrinkage is used. Shrinkage estimation is only available for method ``BSL'' and ``semiBSL''. Current options 
#' are ``glasso'' for the graphical lasso method of Friedman et al (2008) and ``Warton'' for the ridge regularisation method 
#' of Warton (2008).
#' @param penalty		The penalty value to be used for the specified shrinkage method. Must be between zero and one
#' if the shrinkage method is ``Warton''.

#' @param logitTransformBound A \eqn{p} by \eqn{2} numeric matrix indicating the upper and lower bound of parameters if a logit
#' transformation is used on the parameter space, where \eqn{p} is the number of parameters. The default is \code{NULL},
#' which means no logit transformation is used. It is also possible to define other transformations with \code{fnSim}
#' and \code{fnPrior}. The first column contains the lower bound of each parameter and the second column contains the
#' upper bound. Infinite lower or upper bound is also supported, eg. \code{matrix(c(1,Inf,0,10,-Inf,0.5),3,2,byrow=TRUE)}.
#' @param standardise	A logical argument that determines whether to standardise the summary statistics before applying
#' the graphical lasso. This is only valid if method is ``BSL'', shrinkage is ``glasso'' and penalty is not \code{NULL}. The 
#' diagonal elements will not be penalised if the shrinkage method is ``glasso''. The default is \code{FALSE}.
#' @param GRC           A logical argument indicating whether the Gaussian rank correlation matrix (Boudt et al., 2012) 
#' should be used to estimate the covariance matrix in ``BSL'' method. The default is \code{FALSE}, which uses the 
#' sample covariance by default.
#' @param parallel		A logical value indicating whether parallel computing should be used for simulation and summary
#' statistic evaluation. The default is \code{FALSE}. When model simulation is fast, it may be preferable to perform
#' serial computations to avoid significant communication overhead between workers.
#' @param parallelArgs	A list of additional arguments to pass into the \code{foreach} function. Only used when parallel
#' computing is enabled, default is \code{NULL}.
#' @param plotOnTheFly  A logical or numeric argument defining whether or by how many iterations a posterior figure will 
#' be plotted during running. If \code{TRUE}, a plot of approximate univariate posteriors based on the current accepted 
#' samples will be shown by every one thousand iterations. The default is \code{FALSE}.
#' @param verbose       A logical argument indicating whether the iteration numbers (\code{1:M}) and accepted proposal
#' flags should be printed to track progress. The default is \code{FALSE}.
#'
#' @param theta0		Deprecated, will be removed in the future, use \code{model} instead. Initial guess of the parameter 
#' value, which is used as the starting value for MCMC.
#' @param fnSim         Deprecated, will be removed in the future, use \code{model} instead. A function that simulates data 
#' for a given parameter value. The first argument should be the parameters. Other necessary arguments (optional) can be 
#' specified with \code{simArgs}.
#' @param fnSum         Deprecated, will be removed in the future, use \code{model} instead. A function for computing summary 
#' statistics of data. The first argument should be the observed or simulated dataset. Other necessary arguments (optional) 
#' can be specified with \code{sumArgs}.
#' @param fnPrior       Deprecated, will be removed in the future, use \code{model} instead. A function that computes the 
#' log of prior density for a parameter. The default is \code{NULL}, which uses an improper flat prior over the real line 
#' for each parameter. The function must have a single input: a vector of parameter values.
#' @param simArgs	    Deprecated, will be removed in the future, use \code{model} instead. A list of additional arguments 
#' to pass into the simulation function. Only use when the input \code{fnSim} requires additional arguments. The default is 
#' \code{NULL}.
#' @param sumArgs	    Deprecated, will be removed in the future, use \code{model} instead. A list of additional arguments 
#' to pass into the summary statistics function. Only use when the input \code{fnSum} requires additional arguments. The 
# default is \code{NULL}.
#' @param thetaNames	Deprecated, will be removed in the future, use \code{model} instead. A string vector of parameter 
#' names, which must have the same length as the parameter vector. The default is \code{NULL}.
#'
#' @return 				An object of class \code{bsl} is returned, containing the following components:
#' \itemize{
#' \item \code{theta}: MCMC samples from the joint approximate posterior distribution of the parameters.
#' \item \code{loglike}: Accepted MCMC samples of the estimated log-likelihood values.
#' \item \code{call}: The original code that was used to call the method.
#' \item \code{model}: The ``BSLMODEL'' object.
#' \item \code{acceptanceRate}: The acceptance rate of the MCMC algorithm.
#' \item \code{earlyRejectionRate}: The early rejection rate of the algorithm (early rejection may occur when using
#' bounded prior distributions).
#' \item \code{errorRate}: The error rate. If any infinite summary statistic or positive infinite loglike occurs during 
#' the process, it is marked as an error and the proposed parameter will be rejected.
#' \item \code{y}: The input observed data.
#' \item \code{n}: The number of simulations from the model per MCMC iteration.
#' \item \code{M}: The number of MCMC iterations.
#' \item \code{covRandWalk}: The covariance matrix used in multivariate normal random walk proposals.
#' \item \code{method}: The string argument indicating the used method.
#' \item \code{shrinkage}: The string argument indicating the shrinkage method.
#' \item \code{penalty}: The penalty value.
#' \item \code{standardise}: Logical, whether to standardise the summary statistics.
#' \item \code{GRC}: Logical, if Gaussian rank correlation matrix is used.
#' \item \code{logitTransform}: Logical, whether a logit transformation is used in the algorithm.
#' \item \code{logitTransformBound}: The matrix of logitTransformBound.
#' \item \code{parallel}: Logical, whether parallel computing is used in the process.
#' \item \code{parallelArgs}: The list of additional arguments to pass into the \code{foreach} function.
#' \item \code{time}: The running time of class \code{difftime}.
#' }
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
#' Friedman J, Hastie T, & Tibshirani R. (2008). Sparse inverse covariance estimation with the 
#' graphical lasso. Biostatistics.
#' \url{https://doi.org/10.1093/biostatistics/kxm045}
#' 
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://doi.org/10.1198/016214508000000021}
#'
#' Ghurye, S. G., & Olkin, I. (1969).
#' Unbiased estimation of some multivariate probability densities and related functions.
#' The Annals of Mathematical Statistics. \url{https://projecteuclid.org/euclid.aoms/1177697501}
#'
#' Boudt, K., Cornelissen, J., and Croux, C. (2012). The Gaussian rank correlation estimator: 
#' robustness properties. Statistics and Computing, 22(2):471-483.
#'
#' @author    Ziwen An, Leah F. South and Christopher C. Drovandi
#' @seealso   \code{\link{ma2}}, \code{\link{cell}} and \code{\link{mgnk}} for examples. 
#' \code{\link{selectPenalty}} for a function to tune the BSLasso tuning parameter
#' and \code{\link{plot}} for functions related to visualisation.
#' @export
bsl <- function(y, n, M, model, covRandWalk, theta0, fnSim, fnSum, method = c("BSL", "uBSL", 
    "semiBSL")[1], shrinkage = NULL, penalty = NULL, fnPrior = NULL, simArgs = NULL,
    sumArgs = NULL, logitTransformBound = NULL, standardise = FALSE, GRC = FALSE, parallel = FALSE,
    parallelArgs = NULL, thetaNames = NULL, plotOnTheFly = FALSE, verbose = FALSE) {
    if (!method %in% c("BSL", "uBSL", "semiBSL")) {
        stop("method must be one of \"BSL\" or \"uBSL\" or \"semiBSL\"")
    }
    if (!parallel & !is.null(parallelArgs)) {
        warning("\"parallelArgs\" is omitted in serial computing")
    }
    if (is.null(shrinkage) && !is.null(penalty)) {
        warning("\"penalty\" will be ignored since no shrinkage method is specified")
    }
    if (!is.null(shrinkage) && is.null(penalty)) {
        stop("\"penalty\" must be specified to provoke shrinkage method")
    }
    if (standardise && (shrinkage != 'glasso' || method != 'BSL')) {
        warning('standardisation is only supported when method is "BSL" and shrinkage is "glasso"')
    }

	# deprecated arguments
	if (!missing(theta0)) {
	    warning("theta0 will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (!missing(fnSim)) {
	    warning("fnSim will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (!missing(fnSum)) {
	    warning("fnSum will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (!is.null(simArgs)) {
	    warning("simArgs will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (!is.null(sumArgs)) {
	    warning("sumArgs will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (!is.null(fnPrior)) {
	    warning("fnPrior will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (!is.null(thetaNames)) {
	    warning("thetaNames will be deprecated in the future, use model instead, see '?BSLModel'")
	}
	if (missing(model)) {
	    if (is.null(fnPrior)) {
		    fnLogPrior <- NULL
		} else {
		    fnLogPrior <- function(...) log(fnPrior(...))
		}
	    model <- BSLModel(fnSim = fnSim, fnSum = fnSum, simArgs = simArgs, sumArgs = sumArgs,
		             fnLogPrior = fnLogPrior, theta0 = theta0, thetaNames = thetaNames)
	} else {
	    stopifnot(class(model) == 'BSLMODEL')
	}

    p <- length(model@theta0)
	fnLogPrior <- model@fnLogPrior
    logitTransform <- !is.null(logitTransformBound)
    if (logitTransform) {
        if (any(dim(logitTransformBound) != c(p, 2))) {
            stop("\"logitTransformBound\" must be a p by 2 matrix, where p is the length of parameter")
        }
    }

    cl <- match.call()
    startTime <- Sys.time()

	# initialise parameters
    ssy <- do.call(model@fnSum, c(list(y), model@sumArgs))
    ns <- length(ssy)
    thetaCurr <- model@theta0
    loglikeCurr <- Inf
    if (logitTransform) {
        thetaTildeCurr <- paraLogitTransform(thetaCurr, logitTransformBound)
    }
    theta <- array(0, c(M, p), dimnames = list(NULL, thetaNames))
    loglike <- numeric(M)
    countAcc <- countEar <- countErr <- 0

    # map the simulation function
	if (parallel) {
	    myFnSimSum <- function(n, theta) fn(model)$fnPar(n, theta, parallelArgs)
	} else {
	    myFnSimSum <- fn(model)$fn
	}
	
    # plot-on-the-fly
    if (plotOnTheFly) {
	    if (plotOnTheFly == 1) {
		    plotOnTheFly  <- 1000
		}
        oldPar <- par()$mfrow
        a <- floor(sqrt(p))
        b <- ceiling(p/a)
        par(mfrow = c(a, b))
    }

    while (is.infinite(loglikeCurr)) {
        # simulate with thetaProp and calculate the summary statistics
		ssx <- myFnSimSum(n, thetaCurr)
        if (any(is.infinite(ssx))) {
		    stop('Inf detected in the summary statistics vector, this will cause an error in likelihood evaluation')
        }

        # compute the loglikelihood
        loglikeCurr <- switch(method, 
		    BSL = gaussianSynLike(ssy, ssx, shrinkage, penalty, standardise, GRC, log = TRUE, verbose = verbose),
            uBSL = gaussianSynLikeGhuryeOlkin(ssy, ssx, log = TRUE, verbose = verbose), 
			semiBSL = semiparaKernelEstimate(ssy, ssx, shrinkage = shrinkage, penalty = penalty)
			)
    }

    for (i in 1:M) {

        flush.console()
        if (verbose) {
            cat("i =", i, "\n")
        }

        # multivariate normal random walk to the proposed value of theta
        if (!logitTransform) {
            thetaProp <- c(mvtnorm::rmvnorm(1, mean = thetaCurr, sigma = covRandWalk))
            logp2 <- 0
        } else {
            thetaTildeCurr <- paraLogitTransform(thetaCurr, logitTransformBound)
            # thetaTildeProp <- mvrnorm(1, thetaTildeCurr, covRandWalk)
            thetaTildeProp <- mvtnorm::rmvnorm(1, mean = thetaTildeCurr, sigma = covRandWalk)
            thetaProp <- paraLogitBackTransform(thetaTildeProp, logitTransformBound)
            logp2 <- jacobianLogitTransform(thetaTildeProp, logitTransformBound, TRUE) - 
			      jacobianLogitTransform(thetaTildeCurr, logitTransformBound, TRUE)
        }

        # early rejection if the proposed theta falls outside of prior coverage
        # / feasible region
        if (!is.null(fnLogPrior)) {
            logp1 <- fnLogPrior(thetaProp) - fnLogPrior(thetaCurr)
            if (logp1 == -Inf) {
                if (verbose) {
                  cat("*** early rejection ***\n")
                }
                theta[i, ] <- thetaCurr
                loglike[i] <- loglikeCurr
                countEar <- countEar + 1
                next
            }
        } else {
            logp1 <- 0
        }
        prob <- exp(logp1 + logp2)

        # simulate with thetaProp and calculate the summary statistics
		ssx <- myFnSimSum(n, thetaProp)

        # reject if inifite value is detected in ssx
        if (any(is.infinite(ssx))) {
            if (verbose) {
                cat("*** reject (infinite ssx) ***\n")
            }
            theta[i, ] <- thetaCurr
            loglike[i] <- loglikeCurr
            countErr <- countErr + 1
            next
        }

        # compute the loglikelihood
        loglikeProp <- switch(method,
		    BSL = gaussianSynLike(ssy, ssx, shrinkage, penalty, standardise, GRC, log = TRUE, verbose = verbose), 
			uBSL = gaussianSynLikeGhuryeOlkin(ssy, ssx, log = TRUE, verbose = verbose), 
			semiBSL = semiparaKernelEstimate(ssy, ssx, shrinkage = shrinkage, penalty = penalty)
			)

        if (loglikeProp == Inf) {
            if (verbose) {
                cat("*** reject (positive infinite loglike) ***\n")
            }
            theta[i, ] <- thetaCurr
            loglike[i] <- loglikeCurr
            countErr <- countErr + 1
            next
        }

        rloglike <- exp(loglikeProp - loglikeCurr)
        if (runif(1) < prob * rloglike) {
            if (verbose) {
                cat("*** accept ***\n")
            }
            thetaCurr <- thetaProp
            loglikeCurr <- loglikeProp
            countAcc <- countAcc + 1
        }

        theta[i, ] <- thetaCurr
        loglike[i] <- loglikeCurr

        if (plotOnTheFly) {
            if (i %% plotOnTheFly == 0) {
                for (k in 1:p) {
                  plot(density(theta[1:i, k]), main = NA, xlab = thetaNames[k],
                    col = 1, lty = 1)
                }
                # Sys.sleep(0.5)
            }
        }
    }

    accRate <- countAcc/M
    earRate <- countEar/M
    errRate <- countErr/M
    time <- difftime(Sys.time(), startTime)

    if (plotOnTheFly) {
        par(mfrow = oldPar)
    }

	result <- new.bsl(theta = theta, loglike = loglike, call = cl, model = model, 
	    acceptanceRate = accRate, earlyRejectionRate = earRate, errorRate = errRate, 
		y = y, n = n, M = M, covRandWalk = covRandWalk, method = method, 
		shrinkage = shrinkage, penalty = penalty, standardise = standardise, GRC = GRC, 
		logitTransform = logitTransform, logitTransformBound = logitTransformBound, 
		parallel = parallel, parallelArgs = parallelArgs, time = time)
    return(result)
}
