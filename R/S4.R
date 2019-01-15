#' @include S4Functions.R
NULL

setOldClass('difftime')

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("functionOrNULL", c("function", "NULL"))

#' S4 class \code{"bsl"}.
#' @description S4 class \code{"bsl"} for Results of \code{bsl} Function.
#' @rawRd \Rdversion{1.1}
#' @slot theta Object of class \code{"matrix"}. MCMC samples from the joint approximate posterior distribution of the parameters.
#' @slot loglike Object of class \code{"numeric"}. Accepted MCMC samples of the estimated log-likelihood values.
#' @slot acceptanceRate Object of class \code{"numeric"}. The acceptance rate of the MCMC algorithm.
#' @slot earlyRejectionRate Object of class \code{"numeric"}. The early rejection rate of the algorithm (early rejection may occur when using
#' bounded prior distributions).
#' @slot call Object of class \code{"call"}. The original code that was used to call the method.
#' @slot y Object of class \code{"ANY"}. The observed data.
#' @slot n Object of class \code{"numeric"}. The number of simulations from the model per MCMC iteration.
#' @slot M Object of class \code{"numeric"}. The number of MCMC iterations.
#' @slot theta0 Object of class \code{"numeric"}. The initial guess of the parameter value.
#' @slot covRandWalk Object of class \code{"matrix"}. The covariance matrix used in multivariate normal random walk proposals.
#' @slot fnSim Object of class \code{"function"}. The data simulation function.
#' @slot fnSum Object of class \code{"function"}. The function for computing summary statistics of data.
#' @slot method Object of class \code{"character"}. The character argument indicating the used method.
#' @slot shrinkage Object of class \code{"characterOrNULL"}. The character argument indicating the shrinkage method.
#' @slot penalty Object of class \code{"numericOrNULL"}. The penalty value.
#' @slot fnPrior Object of class \code{"functionOrNULL"}. The function that computes the prior density for a parameter.
#' @slot simArgs Object of class \code{"listOrNULL"}. The list of additional arguments to pass into the simulation function.
#' @slot sumArgs Object of class \code{"listOrNULL"}. The list of additional arguments to pass into the summary statistics function.
#' @slot logitTransform Object of class \code{"logical"}. The logical argument indicating whether a logit transformation is used in the algorithm.
#' @slot logitTransformBound Object of class \code{"matrixOrNULL"}. The matrix of logitTransformBound.
#' @slot standardise Object of class \code{"logical"}. The logical argument that determines whether to standardise the summary statistics.
#' @slot parallel Object of class \code{"logical"}. The logical value indicating whether parallel computing is used in the process.
#' @slot parallelArgs Object of class \code{"listOrNULL"}. The list of additional arguments to pass into the \code{foreach} function.
#' @slot thetaNames Object of class \code{"expression"}. The character vector of parameter names.
#' @slot time Object of class \code{"difftime"}. The running time.
#' @export
setClass("bsl", slots = c(theta = "matrix", loglike = "numeric",
                          acceptanceRate = "numeric", earlyRejectionRate = "numeric", call = "call",
                          y = "ANY", n = "numeric", M = "numeric", theta0 = "numeric", covRandWalk = "matrix",
                          fnSim = "function", fnSum = "function", method = "character", shrinkage = "characterOrNULL",
                          penalty = "numericOrNULL", fnPrior = "functionOrNULL", simArgs = "listOrNULL", sumArgs = "listOrNULL",
                          logitTransform = "logical", logitTransformBound = "matrixOrNULL", standardise = "logical",
                          parallel = "logical", parallelArgs = "listOrNULL", thetaNames = "expression",
                          time = "difftime"))

setMethod("initialize", "bsl", initialize.bsl)
setValidity("bsl", check.bsl)

#' Show method for class "bsl". Display the basic information of a bsl object.
#' @param object   A "bsl" class object to be displayed.
#' @exportMethod show
#' @describeIn bsl Display the basic information of a bsl object. See \code{\link{show.bsl}}.
#' @aliases bsl,bsl-method
setMethod('show', signature(object = 'bsl'), show.bsl)

#' Summary method for class "bsl"
#' @description Summarise a bsl class object.
#' @param thetaNames Parameter names to be shown in the summary table. Parameter names of the bsl object will be used by default.
#' @param ... Other arguments.
#' @param y Ignore.
#' @return A vector of the number of simulations per iteration, acceptance rate of the Markov chain annd scaled effective sample size for each parameter.
#' @docType methods
#' @rdname bsl-class
setGeneric('summary')
#' @exportMethod summary
#' @describeIn bsl Summarise a bsl class object. See \code{\link{summary.bsl}}.
#' @aliases bsl,bsl-method
setMethod('summary', 'bsl', summary.bsl)

#' Plot method for class "bsl"
#' @description Plot the univariate marginal posterior plot of a bsl class object.
#' @param x           A "bsl" class object to plot.
#' @param which       An integer argument indicating which plot function to be used. The default, \code{1L}, uses
#' the plain \code{plot} to visualise the result. \code{2L} uses ggplot2 to generate an aesthetically nicer figure.
#' @param thin        A numeric argument indicating the gap between samples to be taken when thinning the MCMC
#' draws. The default is \code{1L}, which means no thinning is used.
#' @param thetaTrue   A set of values to be included on the plots as a reference line. The default is \code{NULL}.
#' @param options.plot  A list of additional arguments to pass into the \code{plot} function. Only use when
#' \code{which} is \code{1L}.
#' @param top         A character argument of the combined plot title if \code{which} is \code{2L}.
#' @param options.density  A list of additional arguments to pass into the \code{geom_density} function. Only use
#' when \code{which} is \code{2L}.
#' @param options.theme  A list of additional arguments to pass into the \code{theme} function. Only use
#' when \code{which} is \code{2L}.
#' @docType methods
#' @rdname bsl-class
setGeneric('plot')
#' @exportMethod plot
#' @describeIn bsl Plot the univariate marginal posterior plot of a bsl class object. See \code{\link{plot.bsl}}.
#' @aliases bsl,bsl-method
setMethod('plot', signature = c(x = 'bsl', y = 'missing'), plot.bsl)
