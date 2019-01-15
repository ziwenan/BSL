#' Performing BSL, BSLasso and semiBSL
#'
#' @description This is the main function for performing MCMC BSL, MCMC BSLasso and MCMC semiBSL.
#' Parallel computing is supported with the R package \code{foreach}.
#'
#' @param y				The observed data - note this should be the raw dataset NOT the set of summary statistics.
#' @param n				The number of simulations from the model per MCMC iteration for estimating the synthetic likelihood.
#' @param M				The number of MCMC iterations.
#' @param theta0		Initial guess of the parameter value, which is used as the starting value for MCMC.
#' @param covRandWalk	A covariance matrix to be used in multivariate normal random walk proposals.
#' @param fnSim         A function that simulates data for a given parameter value. The first argument should be the
#' parameters. Other necessary arguments (optional) can be specified with \code{simArgs}.
#' @param fnSum         A function for computing summary statistics of data. The first argument should be the observed
#' or simulated dataset. Other necessary arguments (optional) can be specified with \code{sumArgs}.
#' @param method        A string argument indicating the method to be used. The default, 'BSL', runs standard BSL or
#' BSLasso if \code{shrinkage} is used. 'semiBSL' runs the semi-parametric BSL algorithm and is more robust to
#' non-normal summary statistics.
#' @param shrinkage      A string argument indicating which shrinkage method to be used. The default is \code{NULL},
#' which means no shrinkage is used. Current options are 'glasso' for graphical lasso and 'Warton' for the
#' ridge regularisation method of Warton (2008).
#' @param penalty		The penalty value to be used for the specified shrinkage method. Must be between zero and one
#' if the shrinkage method is 'Warton'.
#' @param fnPrior		A function that computes the prior density for a parameter. The default is \code{NULL}, which
#' is an improper flat prior over the real line for each parameter. The function must have a single input: a vector
#' of parameter values.
#' @param simArgs	    A list of additional arguments to pass into the simulation function. Only use when the input
#' \code{fnSim} requires additional arguments. The default is \code{NULL}.
#' @param sumArgs	    A list of additional arguments to pass into the summary statistics function. Only use when the
#' input \code{fnSum} requires additional arguments. The default is \code{NULL}.
#' @param logitTransformBound A \eqn{p} by \eqn{2} numeric matrix indicating the upper and lower bound of parameters if a logit
#' transformation is used on the parameter space, where \eqn{p} is the number of parameters. The default is \code{NULL},
#' which means no logit transformation is used. It is also possible to define other transformations with \code{fnSim}
#' and \code{fnPrior}. The first column contains the lower bound of each parameter and the second column contains the
#' upper bound. Infinite lower or upper bound is also supported, eg. \code{matrix(c(1,Inf,0,10,-Inf,0.5),3,2,byrow=TRUE)}.
#' @param standardise	A logical argument that determines whether to standardise the summary statistics before applying
#' the graphical lasso. This is only valid if shrinkage is 'glasso' and penalty is not \code{NULL}. The diagonal
#' elements will not be penalised if the shrinkage method is 'glasso'. The default is \code{FALSE}.
#' @param parallel		A logical value indicating whether parallel computing should be used for simulation and summary
#' statistic evaluation. The default is \code{FALSE}. When model simulation is fast, it may be preferable to perform
#' serial computations to avoid significant communication overhead between workers.
#' @param parallelArgs	A list of additional arguments to pass into the \code{foreach} function. Only used when parallel
#' computing is enabled, default is \code{NULL}.
#' @param thetaNames	A string vector of parameter names, which must have the same length as the parameter vector.
#' The default is \code{NULL}.
#' @param plotOnTheFly  A logical argument. If \code{TRUE}, a plot of approximate univariate posteriors based on the
#' current accepted samples will be shown every 1000 iterations. The default is \code{FALSE}.
#' @param verbose       A logical argument indicating whether the iteration numbers (\code{1:M}) and accepted proposal
#' flags should be printed to track progress. The default is \code{FALSE}.
#'
#' @return 				An object of class \code{bsl} is returned, containing the following components:
#' \itemize{
#' \item \code{theta}: MCMC samples from the joint approximate posterior distribution of the parameters.
#' \item \code{loglike}: Accepted MCMC samples of the estimated log-likelihood values.
#' \item \code{acceptanceRate}: The acceptance rate of the MCMC algorithm.
#' \item \code{earlyRejectionRate}: The early rejection rate of the algorithm (early rejection may occur when using
#' bounded prior distributions).
#' \item \code{call}: The original code that was used to call the method.
#' \item \code{y}: The input observed data.
#' \item \code{n}: The input number of simulations from the model per MCMC iteration.
#' \item \code{M}: The input number of MCMC iterations.
#' \item \code{theta0}: The input initial guess of the parameter value.
#' \item \code{covRandWalk}: The input covariance matrix used in multivariate normal random walk proposals.
#' \item \code{fnSim}: The input data simulation function.
#' \item \code{fnSum}: The input function for computing summary statistics of data.
#' \item \code{method}: The input string argument indicating the used method.
#' \item \code{shrinkage}: The input string argument indicating the shrinkage method.
#' \item \code{penalty}: The input penalty value.
#' \item \code{fnPrior}: The input function that computes the prior density for a parameter.
#' \item \code{simArgs}: The input list of additional arguments to pass into the simulation function.
#' \item \code{sumArgs}: The input list of additional arguments to pass into the summary statistics function.
#' \item \code{logitTransform}: The logical argument indicating whether a logit transformation is used in the algorithm.
#' \item \code{logitTransformBound}: The input matrix of logitTransformBound.
#' \item \code{standardise}: The input logical argument that determines whether to standardise the summary statistics.
#' \item \code{parallel}: The input logical value indicating whether parallel computing is used in the process.
#' \item \code{parallelArgs}: The input list of additional arguments to pass into the \code{foreach} function.
<<<<<<< HEAD
#' \item \code{thetaNames}: The character vector of parameter names.
=======
#' \item \code{thetaNames}: A string vector of parameter names.
>>>>>>> c05a589e2012cf5ea524f46111f3ead1ef61b19c
#' \item \code{time}: The running time of class \code{difftime}.
#' }
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
#' Warton, D. I. (2008). Penalized Normal Likelihood and Ridge Regularization of Correlation and
#' Covariance Matrices, Journal of the American Statistical Association.
#' \url{https://doi.org/10.1198/016214508000000021}
#'
#' @author 								Ziwen An, Leah F. South and Christopher C. Drovandi
#' @seealso 							\code{\link{selectPenalty}} for a function to tune the BSLasso tuning parameter
#' and \code{\link{plot}} for functions related to visualisation.
#' @export
bsl <- function(y, n, M, theta0, covRandWalk, fnSim, fnSum, method = c("BSL",
    "semiBSL")[1], shrinkage = NULL, penalty = NULL, fnPrior = NULL, simArgs = NULL,
    sumArgs = NULL, logitTransformBound = NULL, standardise = FALSE, parallel = FALSE,
    parallelArgs = NULL, thetaNames = NULL, plotOnTheFly = FALSE, verbose = FALSE) {
    if (!method %in% c("BSL", "semiBSL")) {
        stop("method must be either \"BSL\" or \"semiBSL\"")
    }
    if (!parallel & !is.null(parallelArgs)) {
        warning("\"parallelArgs\" is omitted in serial computing")
    }
    if (is.null(shrinkage) && !is.null(penalty)) {
        warning("\"penalty\" will not ignored since no shrinkage method is specified")
    }
    if (!is.null(shrinkage) && is.null(penalty)) {
        stop("\"penalty\" must be specified to provoke shrinkage method")
    }
    if (is.null(penalty) & standardise) {
        warning("standardisation is only supported in BSLasso")
    }

    p <- length(theta0)
    logitTransform <- !is.null(logitTransformBound)
    if (logitTransform) {
        if (any(dim(logitTransformBound) != c(p, 2))) {
            stop("\"logitTransformBound\" must be a p by 2 matrix, where p is the length of parameter")
        }
    }

    cl <- match.call()
    startTime <- Sys.time()
    # match the simulation function
    if (is.null(simArgs)) {
        myFnSim <- function(theta) {
            do.call(fnSim, list(theta))
        }
    } else {
        myFnSim <- function(theta) {
            do.call(fnSim, c(list(theta), simArgs))
        }
    }

    # match the summary statistics function
    if (is.null(sumArgs)) {
        myFnSum <- function(x) {
            do.call(fnSum, list(x))
        }
    } else {
        myFnSum <- function(x) {
            do.call(fnSum, c(list(x), sumArgs))
        }
    }

    if (plotOnTheFly) {
        oldPar <- par()$mfrow
        a <- floor(sqrt(p))
        b <- ceiling(p/a)
        par(mfrow = c(a, b))
    }

    # initialise parameters
    ssy <- myFnSum(y)
    ns <- length(ssy)
    thetaCurr <- theta0
    loglikeCurr <- Inf
    if (logitTransform) {
        thetaTildeCurr <- paraLogitTransform(thetaCurr, logitTransformBound)
    }
    theta <- array(0, c(M, p), dimnames = list(NULL, thetaNames))
    loglike <- numeric(M)
    countAcc <- countEar <- countErr <- 0

    while (is.infinite(loglikeCurr)) {
        # simulate with thetaProp and calculate the summary statistics
        ssx <- Inf
        while (any(is.infinite(ssx))) {
            if (!parallel) {
                ssx <- array(0, c(n, ns))
                for (j in 1:n) {
                    x <- myFnSim(thetaCurr)
                    ssx[j, ] <- myFnSum(x)
                }
            } else {
                # use foreach for parallel computing
                ssx <- do.call(foreach, c(list(j = 1:n, .combine = rbind),
                               parallelArgs)) %dopar% {
                    x <- myFnSim(thetaCurr)
                    myFnSum(x)
                }
            }
        }

        # compute the loglikelihood
        loglikeCurr <- switch(method, BSL = gaussianSynLike(ssy, ssx, shrinkage,
            penalty, standardise, log = TRUE, verbose = verbose), semiBSL = semiparaKernelEstimate(ssy,
            ssx, shrinkage = shrinkage, penalty = penalty))
    }

    for (i in 1:M) {

        flush.console()
        if (verbose) {
            cat("i =", i, "\n")
        }

        # multivariate normal random walk to the proposed value of theta
        if (!logitTransform) {
            thetaProp <- c(rmvnorm(1, mean = thetaCurr, sigma = covRandWalk))
            p2 <- 1
        } else {
            thetaTildeCurr <- paraLogitTransform(thetaCurr, logitTransformBound)
            # thetaTildeProp <- mvrnorm(1, thetaTildeCurr, covRandWalk)
            thetaTildeProp <- rmvnorm(1, mean = thetaTildeCurr, sigma = covRandWalk)
            thetaProp <- paraLogitBackTransform(thetaTildeProp, logitTransformBound)
            p2 <- jacobianLogitTransform(thetaTildeProp, logitTransformBound)/jacobianLogitTransform(thetaTildeCurr,
                logitTransformBound)
        }

        # early rejection if the proposed theta falls outside of prior coverage
        # / feasible region
        if (!is.null(fnPrior)) {
            p1 <- fnPrior(thetaProp)/fnPrior(thetaCurr) # We may want to change this to log prior in a future release for stability
            if (p1 == 0) {
                if (verbose) {
                  cat("*** early rejection ***\n")
                }
                theta[i, ] <- thetaCurr
                loglike[i] <- loglikeCurr
                countEar <- countEar + 1
                next
            }
        } else {
            p1 <- 1
        }
        prob <- p1 * p2

        # simulate with thetaProp and calculate the summary statistics
        if (!parallel) {
            ssx <- array(0, c(n, ns))
            for (j in 1:n) {
                x <- myFnSim(thetaProp)
                ssx[j, ] <- myFnSum(x)
            }
        } else {
            # use foreach for parallel computing
            ssx <- do.call(foreach, c(list(j = 1:n, .combine = rbind),
                parallelArgs)) %dopar% {
                x <- myFnSim(thetaProp)
                myFnSum(x)
            }
        }

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
        loglikeProp <- switch(method, BSL = gaussianSynLike(ssy, ssx, shrinkage,
            penalty, standardise, log = TRUE, verbose = verbose), semiBSL = semiparaKernelEstimate(ssy,
            ssx, shrinkage = shrinkage, penalty = penalty))

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

        # cat('theta: ',thetaProp,'\t') cat('ll: ',loglikeProp,'\t') cat('p:
        # ',p,'\t') cat('prob: ',p * rloglike,'\n') accept the proposed theta
        # with a probability
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
            if (i %% 1000 == 0) {
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

    results <- new('bsl', theta = theta, loglike = loglike,
        acceptanceRate = accRate, earlyRejectionRate = earRate,
        call = cl, y = y, n = n, M = M, theta0 = theta0, covRandWalk = covRandWalk,
        fnSim = fnSim, fnSum = fnSum, method = method, shrinkage = shrinkage,
        penalty = penalty, fnPrior = fnPrior, simArgs = simArgs, sumArgs = sumArgs,
        logitTransform = logitTransform, logitTransformBound = logitTransformBound,
        standardise = standardise, parallel = parallel, parallelArgs = parallelArgs,
        thetaNames = as.expression(thetaNames), time = time)
    return(results)
}
