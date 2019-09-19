setClassUnion("functionOrNULL", c("function", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))

#' S4 class ``BSLMODEL''
#' @description The S4 class contains the simulation and summary statistics function and other necessary arguments
#' for a model to run in the main \code{bsl} function.
#' @rawRd \Rdversion{1.1}
#' @slot fnSim A function that simulates data for a given parameter value. The first argument should be the
#' parameters. Other necessary arguments (optional) can be specified with \code{simArgs}.
#' @slot fnSimVec A vectorised function that simulates a number of datasets simultaneously for a given parameter
#' value. If this is not \code{NULL}, vectorised simulation function will be used instead of \code{fnSim}. The
#' first two arguments should be the number of simulations to run and parameters, respectively. Other necessary 
#' arguments (optional) can be specified with \code{simArgs}. The output must be a list of each simulation result.
#' @slot fnSum A function for computing summary statistics of data. The first argument should be the observed
#' or simulated dataset. Other necessary arguments (optional) can be specified with \code{sumArgs}. The users
#' should code this function carefully so the output have fixed length and never contain any \code{Inf} value.
#' @slot fnLogPrior A function that computes the log of prior density for a parameter. The default is \code{NULL},
#' which uses an improper flat prior over the real line for each parameter. The function must have a single input: a
#' vector of parameter values.
#' @slot simArgs A list of additional arguments to pass into the simulation function. Only use when the input
#' \code{fnSim} or \code{fnSimVec} requires additional arguments. The default is \code{NULL}.
#' @slot sumArgs A list of additional arguments to pass into the summary statistics function. Only use when the
#' input \code{fnSum} requires additional arguments. The default is \code{NULL}.
#' @slot theta0 Initial guess of the parameter value, which is used as the starting value for MCMC.
#' @slot thetaNames Expression, parameter names.
#' @slot ns The number of summary statistics of a single observation. Note this will be generated automatically, 
#' thus is not required for initialisation.
#' @slot test Logical indicator of whether a short simulation test will be ran upon initialisation.
#' @exportClass BSLMODEL
setClass('BSLMODEL', slots = c(fnSim = 'functionOrNULL', fnSimVec = 'functionOrNULL', fnSum = 'functionOrNULL',
                                           fnLogPrior = 'functionOrNULL', simArgs = 'listOrNULL', sumArgs = 'listOrNULL',
                                           theta0 = 'numeric', thetaNames = 'expression', ns = 'integer', test = 'logical'))

#' Constructor for class ``BSLMODEL''
#' @description \code{BSLModel} is the constructor function for a \code{BSLMODEL} object.
#' @examples
#' data(ma2)
#' model <- BSLModel(fnSim = ma2_sim, fnSum = ma2_sum, simArgs = ma2$sim_options, theta0 = ma2$start,
#'                   fnLogPrior = ma2_logPrior)
#' validObject(model)
#' @rdname BSLMODEL-class
#' @export BSLModel
BSLModel <- function(fnSim, fnSimVec, fnSum, fnLogPrior, simArgs, sumArgs, theta0, thetaNames, test = TRUE) {
    new(Class = "BSLMODEL", fnSim = fnSim, fnSimVec = fnSimVec, fnSum = fnSum, fnLogPrior = fnLogPrior, 
	    simArgs = simArgs, sumArgs = sumArgs, theta0 = theta0, thetaNames = thetaNames, test = test)
}

#' @description \code{initialize.BSLMODEL} initialises a ``BSLMODEL'' object.
#' @param .Object       A ``BSLMODEL'' object.
#' @param fnSim         A function that simulates data for a given parameter value. The first argument should be the 
#' parameters. Other necessary arguments (optional) can be specified with \code{simArgs}.
#' @param fnSimVec      A vectorised function that simulates a number of datasets simultaneously for a given parameter
#' value. The first two arguments should be the number of simulations to run and parameters, respectively. Other necessary 
#' arguments (optional) can be specified with \code{simArgs}. The output must be a list of each simulation result.
#' @param simArgs	    A list of additional arguments to pass into the simulation function. Only use when the input 
#' \code{fnSim} requires additional arguments.
#' @param fnSum         A function for computing summary statistics of data. The first argument should be the observed 
#' or simulated dataset. Other necessary arguments (optional) can be specified with \code{sumArgs}.
#' @param sumArgs	    A list of additional arguments to pass into the summary statistics function. Only use when the 
#' input \code{fnSum} requires additional arguments.
#' @param fnLogPrior    A function that computes the log of prior density for a parameter. If this is missing, the prior 
#' by default is an improper flat prior over the real line for each parameter. The function must have a single input: a
#' vector of parameter values.
#' @param theta0        Initial guess of the parameter value.
#' @param thetaNames	A string vector of parameter names, which must have the same length as the parameter vector. 
#' @param test          Logical, whether a short simulation test will be ran upon initialisation.
#' @rdname BSLMODEL-class
#' @aliases BSLModel
initialize.BSLMODEL <- function(.Object, fnSim, fnSimVec, simArgs, fnSum, sumArgs, fnLogPrior, theta0, thetaNames, 
    test = TRUE) {
    cat('*** initialize "BSLMODEL" ***\n')
    has.fnSim <- !missing(fnSim) || !missing(fnSimVec)
    cat(paste('has simulation function:', has.fnSim, '\n'))
    has.fnSum <- !missing(fnSum)
    cat(paste('has summary statistics function:', has.fnSum, '\n'))
    has.theta0 <- !missing(theta0)
    cat(paste('has initial guess / point estimate of the parameter:', has.theta0, '\n'))
    if (has.fnSim && has.fnSum && has.theta0) {
        if (!missing(fnSim)) .Object@fnSim <- fnSim
        if (!missing(fnSimVec)) .Object@fnSimVec <- fnSimVec
        if (!missing(simArgs)) .Object@simArgs <- simArgs
        .Object@fnSum <- fnSum
        if (!missing(sumArgs)) .Object@sumArgs <- sumArgs
        if (!missing(fnLogPrior)) .Object@fnLogPrior <- fnLogPrior
        if (is.null(.Object@fnLogPrior)) {
            .Object@fnLogPrior <- function(theta) 0
            cat('No prior has been defined in the model, use the default improper flat prior\n')
        }
        .Object@theta0 <- theta0
		.Object@test <- test
        validObject(.Object)
        .Object@ns <- setns(.Object)

        if (missing(thetaNames)) { # missing
            if (!is.null(names(theta0))) { # use the name of theta0
                .Object@thetaNames <- as.expression(names(theta0))
            } else { # use default theta names
                thetaNames <- vector('expression', length(theta0))
                for (i in 1 : length(theta0)) {
                    thetaNames[i] <- as.expression(substitute(theta[j], list(j = i)))
                }
                .Object@thetaNames <- thetaNames
            }
        } else { # !missing
            if (is.null(thetaNames)) { # use default theta names
                thetaNames <- vector('expression', length(theta0))
                for (i in 1 : length(theta0)) {
                    thetaNames[i] <- as.expression(substitute(theta[j], list(j = i)))
                }
                .Object@thetaNames <- thetaNames
            } else { # has thetaNames
                if (length(thetaNames) != length(theta0)) {
                    cat(paste('The length of thetaNames does not match the length of theta0,', length(theta0), '\n'))
                } else {
                    .Object@thetaNames <- as.expression(thetaNames)
                }
            }
        }
    } else {
        cat('an empty (invalid) BSLMODEL object has been created due to one or more missing slots\n')
    }
    cat('*** end initialize ***\n')
    return (.Object)
}

#' @description \code{validBSLModelObject} check the validity of a ``BSLMODEL'' object.
#' @param object       A ``BSLMODEL'' object.
#' @rdname BSLMODEL-class
#' @aliases BSLModel
validBSLModelObject <- function(object) {
    if (is.null(object@fnSim) && is.null(object@fnSimVec)) {
        return('No available simulation function is provided')
    }
	
	if (object@test) {
	    cat('running a short simulation test ... \n')
	    # test simulation function
        if (!is.null(object@fnSimVec)) {
            x <- try(do.call(object@fnSimVec, c(list(10, object@theta0), object@simArgs)))
            if (inherits(x, 'try-error')) {
                return('Fail to run simulations with the given vectorised simulation function')
            }
            if (!is.list(x) || length(x) != 10) {
                return('The output of the given vectorised simulation function must be a list of the
                       same length as the number of iterations')
            }
            } else {
                x <- list()
                x[[1]] <- try(do.call(object@fnSim, c(list(object@theta0), object@simArgs)))
                if (inherits(x[[1]], 'try-error')) {
                    return('Fail to run simulations with the given simulation function')
                }
        }
        
        # test summary statistics function
        ssx <- try(do.call(object@fnSum, c(list(x[[1]]), object@sumArgs)))
        if (inherits(ssx, 'try-error')) {
            return('Fail to get summary statistics with the given summary statistics function')
        }
        if (!is.vector(ssx)) {
            return('The output of the summary statistics function must be a vector')
        }
	} else {
	    cat('omitting short simulation test\n')
	}
    
    # prior
    if (object@fnLogPrior(object@theta0) == -Inf) {
        return('The given parameter value theta0 has no prior support\n')
    }

    # pass all checks
    TRUE
}

setMethod("initialize", "BSLMODEL", initialize.BSLMODEL)
setValidity("BSLMODEL", validBSLModelObject)

fnBSLModel <- function(.Object) {
    if (!is.null(.Object@fnSimVec)) { # use vectorised simulation function
        fnPar <- function(n, theta, parallelArgs = list()) {
		    j <- NULL
            parallelArgs$.export <- c(parallelArgs$.export, '.Object')
            x <- do.call(.Object@fnSimVec, c(list(n, theta), .Object@simArgs))
            do.call(foreach, c(list(j = 1:n, .combine = rbind), parallelArgs)) %dopar% {
                do.call(.Object@fnSum, c(list(x[[j]]), .Object@sumArgs))
            }
        }
        fn <- function(n, theta) {
            x <- do.call(.Object@fnSimVec, c(list(n, theta), .Object@simArgs))
            matrix(sapply(x, FUN = function(y) do.call(.Object@fnSum, c(list(y), .Object@sumArgs))),
                   nrow = n, byrow = TRUE)
        }
    } else { # non-vectorised simulation function
        fnPar <- function(n, theta, parallelArgs = list()) {
		    j <- NULL
            parallelArgs$.export <- c(parallelArgs$.export, '.Object')
            do.call(foreach, c(list(j = 1:n, .combine = rbind), parallelArgs)) %dopar% {
                x <- do.call(.Object@fnSim, c(list(theta), .Object@simArgs))
                do.call(.Object@fnSum, c(list(x), .Object@sumArgs))
            }
        }

        x <- do.call(.Object@fnSim, c(list(.Object@theta0), .Object@simArgs))
        ns <- length(do.call(.Object@fnSum, c(list(x), .Object@sumArgs)))
        fn <- function(n, theta) {
            ssx <- array(0, c(n, ns))
            for (j in 1:n) {
                x <- do.call(.Object@fnSim, c(list(theta), .Object@simArgs))
                ssx[j, ] <- do.call(.Object@fnSum, c(list(x), .Object@sumArgs))
            }
            ssx
        }
    }
    return(list(fn = fn, fnPar = fnPar))
}

#' Functions to be used in bsl (for internal use)
#' @description Generate the generic simulation and summary statistics function for \code{n} simulations
#' and fixed \code{theta} (for internal use).
#' @inheritParams validBSLModelObject
setGeneric('fn', function(.Object) standardGeneric("fn"))
#' @describeIn BSLMODEL Generate the generic simulation and summary statistics function for \code{n}
#' simulations and fixed \code{theta} (for internal use).
#' @exportMethod fn
setMethod('fn', signature = c(.Object = 'BSLMODEL'), fnBSLModel)



setLengthSummStat <- function(.Object) {
    if (!is.null(.Object@fnSimVec)) {
        x <- do.call(.Object@fnSimVec, c(list(2, .Object@theta0), .Object@simArgs))
    } else {
        x <- list()
        x[[1]] <- do.call(.Object@fnSim, c(list(.Object@theta0), .Object@simArgs))
    }
    ns <- length(do.call(.Object@fnSum, c(list(x[[1]]), .Object@sumArgs)))
    return(as.integer(ns))
}

#' Find the length of summary statistics (for internal use)
#' @description Find and set the length of summary statistics with a test run (for internal use).
#' @inheritParams validBSLModelObject
setGeneric('setns', valueClass = "integer", function(.Object) standardGeneric('setns'))
#' @describeIn BSLMODEL Find and set the length of summary statistics with a test run (for internal use).
#' @exportMethod setns
setMethod('setns', signature = c(.Object = 'BSLMODEL'), setLengthSummStat)
