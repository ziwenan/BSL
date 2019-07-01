#' @include S4Functions.R
NULL

#' @include BSL_model.R
NULL

setOldClass('difftime')

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
# setClassUnion("listOrNULL", c("list", "NULL")) # defined in BSL_model.R
# setClassUnion("functionOrNULL", c("function", "NULL")) # defined in BSL_model.R

#' S4 class ``bsl''.
#' @description The result from function \code{bsl} is saved as class ``BSL''.
#' @rawRd \Rdversion{1.1}
#' @slot theta Object of class ``matrix''. MCMC samples from the joint approximate posterior distribution of the parameters.
#' @slot loglike Object of class ``numeric''. Accepted MCMC samples of the estimated log-likelihood values.
#' @slot call Object of class ``call''. The original code that was used to call the method.
#' @slot model Object of class ``BSLMODEL''.
#' @slot acceptanceRate Object of class ``numeric''. The acceptance rate of the MCMC algorithm.
#' @slot earlyRejectionRate Object of class ``numeric''. The early rejection rate of the algorithm (early rejection may occur when using
#' bounded prior distributions).
#' @slot errorRate Object of class ``numeric''. The error rate. If any infinite summary statistic or positive infinite loglike occurs during
#' the process, it is marked as an error and the proposed parameter will be rejected.
#' @slot y Object of class ``ANY''. The observed data.
#' @slot n Object of class ``numeric''. The number of simulations from the model per MCMC iteration.
#' @slot M Object of class ``numeric''. The number of MCMC iterations.
#' @slot covRandWalk Object of class ``matrix''. The covariance matrix used in multivariate normal random walk proposals.
#' @slot method Object of class ``character''. The character argument indicating the used method.
#' @slot shrinkage Object of class ``characterOrNULL''. The character argument indicating the shrinkage method.
#' @slot penalty Object of class ``numericOrNULL''. The penalty value.
#' @slot GRC Object of class ``logical''. Whether the Gaussian rank correlation matrix is used.
#' @slot logitTransform Object of class ``logical''. The logical argument indicating whether a logit transformation is used in the algorithm.
#' @slot logitTransformBound Object of class ``matrixOrNULL''. The matrix of logitTransformBound.
#' @slot standardise Object of class ``logical''. The logical argument that determines whether to standardise the summary statistics.
#' @slot parallel Object of class ``logical''. The logical value indicating whether parallel computing is used in the process.
#' @slot parallelArgs Object of class ``listOrNULL''. The list of additional arguments to pass into the \code{foreach} function.
#' @slot time Object of class ``difftime''. The running time.
#' @export new.bsl
#' @exportClass bsl
new.bsl <- setClass("bsl", slots = c(theta = "matrix", loglike = "numeric", call = "call", model = 'BSLMODEL', 
                          acceptanceRate = "numeric", earlyRejectionRate = "numeric", errorRate = "numeric", 
                          y = "ANY", n = "numeric", M = "numeric", covRandWalk = "matrix",
                          method = "character", shrinkage = "characterOrNULL", penalty = "numericOrNULL", 
						  standardise = "logical", GRC = 'logical', 
                          logitTransform = "logical", logitTransformBound = "matrixOrNULL", 
                          parallel = "logical", parallelArgs = "listOrNULL", time = "difftime"))

setValidity("bsl", check.bsl)

#' Show method for class ``bsl''. Display the basic information of a bsl object.
#' @param object   A ``bsl'' class object to be displayed.
#' @exportMethod show
#' @describeIn bsl Display the basic information of a ``bsl'' object. See \code{\link{show.bsl}}.
#' @aliases bsl,bsl-method
setMethod('show', signature(object = 'bsl'), show.bsl)

#' Summary method for class ``bsl''
#' @description Summarise a ``bsl'' class object.
#' @param thetaNames Parameter names to be shown in the summary table. If not given, parameter names of the ``bsl'' object will be used by default.
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

#' Plot method for class ``bsl''
#' @description Plot the univariate marginal posterior plot of a ``bsl'' class object.
#' @param x           A ``bsl'' class object to plot.
#' @param which       An integer argument indicating which plot function to be used. The default, \code{1L}, uses
#' the plain \code{plot} to visualise the result. \code{2L} uses ggplot2 to draw the plot.
#' @param thin        A numeric argument indicating the gap between samples to be taken when thinning the MCMC
#' draws. The default is \code{1}, which means no thinning is used.
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
#' @describeIn bsl Plot the univariate marginal posterior plot of a ``bsl'' class object. See \code{\link{plot.bsl}}.
#' @aliases bsl,bsl-method
setMethod('plot', signature = c(x = 'bsl', y = 'missing'), plot.bsl)
