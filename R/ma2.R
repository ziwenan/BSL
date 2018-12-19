#' An MA(2) model
#'
#' @description In this example we wish to estimate the parameters of a simple MA(2) time series model.
#' We provide the data and tuning parameters required to reproduce the results in An et al. (2018).
#'
#' @param theta     A vector of proposed model parameters, \eqn{\theta_1} and \eqn{\theta_2}
#' @param options	A list of options for simulating data from the model. For this example, the list contains only \eqn{T}, the number of observations.
#' @param x			Observed or simulated data in the format of a vector of length \eqn{T}.
#' @param T         The length of the MA(2) to be simulated.
#'
#' @details
#' This example is based on estimating the parameters of a basic MA(2) time series model
#' of the form
#'
#' \deqn{y_t = z_t + \theta_1 z_{t-1} + \theta_2 z_{t-2},}
#'

#' where \eqn{t=1,\ldots,T} and \eqn{z_t ~ N(0,1)} for \eqn{t=-1,0,\ldots,T}.
#' A uniform prior is used for this example, subject to the restrictions that
#' \eqn{-2<\theta_1<2}, \eqn{\theta_1+\theta_2>-1} and \eqn{\theta_1-\theta_2<1}
#' so that invertibility of the time series is satisfied. The summary statistics
#' are simply the full data.
#'
#' @section A simulated dataset:
#'
#' An example 'observed' dataset and the tuning parameters relevant to that example
#' can be obtained using \code{data(ma2)}. This 'observed' data is a simulated dataset
#' with \eqn{\theta_1 = 0.6}, \eqn{\theta_2=0.2} and \eqn{T=50}. Further information
#' about this model and the specific choices of tuning parameters used in BSL and
#' BSLasso can be found in An et al. (2018).
#'
#' \itemize{
#'  \item \code{data}: A time series dataset, in the form of a vector of length \eqn{T}
#'  \item \code{sim_options}: A list containing \eqn{T=50}
#'  \item \code{start}: A vector of suitable initial values of the parameters for MCMC
#'  \item \code{cov}: Covariance matrix of the multivariate normal random walk, in the form of a \eqn{2 x 2} matrix
#' }
#'
#' @examples
#' # Loading the data for this example
#' data(ma2)
#' true_ma2 <- c(0.6,0.2)
#'
#' # Performing BSL (reduce the number of iterations M if desired)
#' resultMa2BSL <- bsl(y = ma2$data, n = 500, M = 300000, theta0 = ma2$start, covRandWalk = ma2$cov,
#'                     fnSim = ma2_sim, fnSum = ma2_sum, fnPrior = ma2_prior, simArgs = ma2$sim_options,
#'                     thetaNames = expression(theta[1], theta[2]), verbose = TRUE)
#' show(resultMa2BSL)
#' summary(resultMa2BSL)
#' plot(resultMa2BSL, which = 1, thetaTrue = true_ma2, thin = 20)
#'
#' # Performing tuning for BSLasso
#' lambda_all <- list(exp(seq(-3,0.5,length.out=20)), exp(seq(-4,-0.5,length.out=20)),
#'                    exp(seq(-5.5,-1.5,length.out=20)), exp(seq(-7,-2,length.out=20)))
#' sp_ma2 <- selectPenalty(ssy = ma2_sum(ma2$data), n = c(50, 150, 300, 500), lambda_all,
#'                         theta = true_ma2, M = 100, sigma = 1.5, fnSim = ma2_sim,
#'                         fnSum = ma2_sum, simArgs = ma2$sim_options)
#' sp_ma2
#' plot(sp_ma2)
#'
#' # Performing BSLasso with a fixed penalty (reduce the number of iterations M if desired)
#' resultMa2BSLasso <- bsl(y = ma2$data, n = 300, M = 250000, theta0 = ma2$start, covRandWalk = ma2$cov,
#'                         fnSim = ma2_sim, fnSum = ma2_sum, shrinkage = 'glasso', penalty = 0.027,
#'                         fnPrior = ma2_prior, simArgs = ma2$sim_options,
#'                         thetaNames = expression(theta[1], theta[2]), verbose = TRUE)
#' show(resultMa2BSLasso)
#' summary(resultMa2BSLasso)
#' plot(resultMa2BSLasso, which = 1, thetaTrue = true_ma2, thin = 20)
#'
#' # Performing semiBSL (reduce the number of iterations M if desired)
#' resultMa2SemiBSL <- bsl(y = ma2$data, n = 500, M = 300000, theta0 = ma2$start, covRandWalk = ma2$cov,
#'                         fnSim = ma2_sim, fnSum = ma2_sum, method = 'semiBSL', fnPrior = ma2_prior,
#'                         simArgs = ma2$sim_options, thetaNames = expression(theta[1], theta[2]),
#'                         verbose = TRUE)
#' show(resultMa2SemiBSL)
#' summary(resultMa2SemiBSL)
#' plot(resultMa2SemiBSL, which = 1, thetaTrue = true_ma2, thin = 20)
#'
#' # Plotting the results together for comparison
#' # plot using the R default plot function
#' par(mar = c(5, 4, 1, 2), oma = c(0, 1, 2, 0))
#' combinePlotsBSL(list(resultMa2BSL, resultMa2BSLasso, resultMa2SemiBSL), which = 1, thetaTrue = true_ma2, thin = 20,
#'                 label = c('bsl', 'bslasso', 'semiBSL'), col = c('red', 'blue', 'green'), lty = 2:4, lwd = 1)
#' mtext('Approximate Univariate Posteriors', outer = TRUE, cex = 1.5)
#'
#' # plot using the ggplot2 package
#' combinePlotsBSL(list(resultMa2BSL, resultMa2BSLasso, resultMa2SemiBSL), which = 2, thetaTrue = true_ma2, thin = 20,
#'                 label = c('bsl   ', 'bslasso   ', 'semiBSL'), options.color = list(values=c('red', 'blue', 'green')),
#'                 options.linetype = list(values = 2:4), options.size = list(values = rep(1, 3)),
#'                 options.theme = list(plot.margin = unit(rep(0.03,4), "npc"), axis.title = element_text(size = 12),
#'                                      axis.text = element_text(size = 8), legend.text = element_text(size = 12)))
#'
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating Bayesian synthetic
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#'
#' @author Ziwen An, Christopher C. Drovandi and Leah F. South
#'
#' @name ma2
#' @usage data(ma2)
NULL


#' The function \code{ma2_sim(theta,T)} simulates an MA(2) time series.
#' @rdname ma2
#' @export
ma2_sim <- function(theta, T) {
    rand <- rnorm(T + 2)
    y <- rand[3 : (T+2)] + theta[1] * rand[2 : (T+1)] + theta[2] * rand[1 : T]
    return(y)
}

#' The function \code{ma2_sum(x)} returns the summary statistics for a given data set. Since the summary statistics are the data, this function simply returns the data.
#' @rdname ma2
#' @export
ma2_sum <- function(x) {
    return(x)
}

#' \code{ma2_prior(theta)} evaluates the (unnormalised) prior, which is uniform subject to several restrictions related to invertibility of the time series.
#' @rdname ma2
#' @export
ma2_prior <- function(theta) {
    theta[2] < 1 & sum(theta) > -1 & diff(theta) > -1
}
