#' Cell biology example
#'
#' @description This example estimates the probabilities of cell motility and
#' cell proliferation for a discrete-time stochastic model of
#' cell spreading. We provide the data and tuning parameters required to
#' reproduce the results in An et al. (2018).
#'
#' @param theta	   A vector of proposed model parameters, \eqn{Pm} and \eqn{Pp}.
#' @param Y        A \code{rows} \eqn{x} \code{cols} \eqn{x} \code{num_obs} array of the cell presences at times \code{1:num_obs} (not time 0).
#' @param Yinit    The initial matrix of cell presences of size \code{rows} \eqn{x} \code{cols}.
#' @param rows     The number of rows in the lattice (rows in the cell location matrix).
#' @param cols     The number of columns in the lattice (columns in the cell location matrix).
#' @param sim_iters The number of discretisation steps to get to when an observation is
#' actually taken. For example, if observations are taken every 5 minutes but the discretisation level is 2.5 minutes,
#' then \code{sim_iters} would be 2. Larger values of \code{sim_iters} lead to more "accurate" simulations from the model, but they also increase the simulation time.
#' @param num_obs  The total number of images taken after initialisation.
#
#' @details
#' Cell motility (movement) and proliferation (reproduction) cause
#' tumors to spread and wounds to heal. If we can measure cell proliferation
#' and cell motility under different situations, then we may be able to use
#' this information to determine the efficacy of different medical treatments.
#'
#' A common method for measuring in vitro cell movement and proliferation is
#' the scratch assay. Cells form a layer on an assay and, once
#' they are completely covering the assay, a scratch is
#' made to separate the cells. Images of the cells are taken until the
#' scratch has closed up and the cells are in contact again.
#' Each image can be converted to a binary matrix by forming a lattice
#' and recording the binary matrix (of size \code{rows} \eqn{x} \code{cols}) of cell presences.
#'
#' The model that we consider is a random walk model with parameters for the probability
#' of cell movement (\eqn{Pm}) and the probability of cell proliferation (\code{Pp})
#' and it has no tractable likelihood function. We use the uninformative priors
#' \eqn{\theta_1 ~ U(0,1)} and \eqn{\theta_2 ~ U(0,1)}.
#'
#' We have a total of 145 summary statistics, which are made up of the Hamming distances
#' between the binary matrices for each time point and the total number of cells at the final time.
#'
#' Details about the types of cells that this model is suitable for
#' and other information can be found in Price et al. (2018) and An et al. (2018). Johnston et al. (2014)
#' use a different ABC method and different summary statistics for a similar example.
#'
#' @section A simulated dataset:
#'
#' An example 'observed' dataset and the tuning parameters relevant to that example
#' can be obtained using \code{data(cell)}. This 'observed' data is a simulated dataset
#' with \eqn{Pm = 0.35} and \eqn{Pp = 0.001}. The lattice has 27 \code{rows} and 36 \code{cols}
#' and there are \code{num_obs = 144} observations after time 0
#' (to mimic images being taken every 5 minutes for 12 hours).
#' The simulation is based on there initially being 110 cells in the assay.
#'
#' Further information about the specific choices of tuning parameters
#' used in BSL and BSLasso can be found in An et al. (2018).
#'
#' \itemize{
#'  \item \code{data}:  The \code{rows} \eqn{x} \code{cols} \eqn{x} \code{num_obs} array of the cell presences at times 1:144.
#'  \item \code{sim_options}: Values of \code{sim_options} relevant to this example.
#'  \item \code{sum_options}: Values of \code{sim_options} relevant to this example, i.e. just the value of \code{Yinit}.
#'  \item \code{start}: A vector of suitable initial values of the parameters for MCMC.
#'  \item \code{cov}: Covariance matrix of the multivariate normal random walk, in the form of a \eqn{2 x 2} matrix.
#' }
#'
#' @examples
#' \donttest{
#' require(doParallel) # You can use a different package to set up the parallel backend
#'
#' # Loading the data for this example
#' data(cell)
#' true_cell <- c(0.35, 0.001)
#'
#' # Performing BSL (reduce the number of iterations M if desired)
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' resultCellBSL <- bsl(cell$data, n = 5000, M = 10000, theta0 = cell$start, covRandWalk = cell$cov,
#'                      fnSim = cell_sim, fnSum = cell_sum, fnPrior = cell_prior,
#'                      simArgs = cell$sim_options, sumArgs = cell$sum_options,
#'                      parallel = TRUE, parallelArgs = list(.packages = 'BSL'),
#'                      thetaNames = expression(P[m], P[p]), verbose = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' show(resultCellBSL)
#' summary(resultCellBSL)
#' plot(resultCellBSL, thetaTrue = true_cell, thin = 20)
#'
#' # Performing tuning for BSLasso
#' lambda_all <- list(exp(seq(0.5,2.5,length.out=20)), exp(seq(0,2,length.out=20)),
#'                    exp(seq(-1,1,length.out=20)), exp(seq(-1,1,length.out=20)))
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' sp_cell <- selectPenalty(ssy = cell_sum(cell$data, cell$sum_options$Yinit),
#'                          n = c(500, 1000, 1500, 2000), lambda_all, theta = true_cell,
#'                          M = 100, sigma = 1.5, fnSim = cell_sim,
#'                          fnSum = cell_sum, simArgs = cell$sim_options,
#'                          sumArgs = cell$sum_options, parallelSim = TRUE,
#'                          parallelSimArgs = list(.packages = 'BSL'), parallelMain = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' sp_cell
#' plot(sp_cell)
#'
#' # Performing BSLasso with a fixed penalty (reduce the number of iterations M if desired)
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' resultCellBSLasso <- bsl(cell$data, n = 1500, M = 10, theta0 = cell$start,
#'                          covRandWalk = cell$cov, fnSim = cell_sim, fnSum = cell_sum,
#'                          shrinkage = 'glasso', penalty = 1.3, fnPrior = cell_prior,
#'                          simArgs = cell$sim_options, sumArgs = cell$sum_options,
#'                          parallel = TRUE, parallelArgs = list(.packages = 'BSL'),
#'                          thetaNames = expression(P[m], P[p]), verbose = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' show(resultCellBSLasso)
#' summary(resultCellBSLasso)
#' plot(resultCellBSLasso, thetaTrue = true_cell, thin = 20)
#'
#' # Performing semiBSL (reduce the number of iterations M if desired)
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' resultCellSemiBSL <- bsl(cell$data, n = 5000, M = 10000, theta0 = cell$start,
#'                          covRandWalk = cell$cov, fnSim = cell_sim, fnSum = cell_sum,
#'                          method = 'semiBSL', fnPrior = cell_prior,
#'                          simArgs = cell$sim_options, sumArgs = cell$sum_options,
#'                          parallel = TRUE, parallelArgs = list(.packages = 'BSL'),
#'                          thetaNames = expression(P[m], P[p]), verbose = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' show(resultCellSemiBSL)
#' summary(resultCellSemiBSL)
#' plot(resultCellSemiBSL, thetaTrue = true_cell, thin = 20)
#'
#' # Plotting the results together for comparison
#' # plot using the R default plot function
#' par(mar = c(5, 4, 1, 2), oma = c(0, 1, 2, 0))
#' combinePlotsBSL(list(resultCellBSL, resultCellBSLasso, resultCellSemiBSL),
#'     which = 1, thetaTrue = true_cell, thin = 20, label = c('bsl', 'bslasso', 'semiBSL'),
#'     col = c('red', 'blue', 'green'), lty = 2:4, lwd = 1)
#' mtext('Approximate Univariate Posteriors', outer = TRUE, cex = 1.5)
#'
#' # plot using the ggplot2 package
#' combinePlotsBSL(list(resultCellBSL, resultCellBSLasso, resultCellSemiBSL),
#'     which = 2, thetaTrue = true_cell, thin = 20, label = c('bsl   ', 'bslasso   ', 'semiBSL'),
#'     options.color = list(values=c('red', 'blue', 'green')),
#'     options.linetype = list(values = 2:4), options.size = list(values = rep(1, 3)),
#'     options.theme = list(plot.margin = grid::unit(rep(0.03,4),"npc"),
#'     axis.title = ggplot2::element_text(size=12), axis.text = ggplot2::element_text(size = 8), 
#'     legend.text = ggplot2::element_text(size = 12)))
#'
#' }
#'
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating Bayesian synthetic
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#'
#' Johnston, S., Simpson, M. J., McElwain, D. L. S., Binder, B. J. &
#' Ross, J. V. (2014). Interpreting Scratch Assays Using Pair Density
#' Dynamic and Approximate Bayesian Computation. Open Biology, 4, 1-11.
#'
#' Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018).
#' Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2017.1302882}
#'
#' @author 								Ziwen An, Leah F. South and Christopher C. Drovandi
#' @name cell
NULL



#' The function \code{cell_sim(theta, Yinit, rows, cols, sim_iters, num_obs)} simulates data from the model, using C++ in the backend.
#' @rdname cell
#' @export
cell_sim <-function(theta, Yinit, rows, cols, sim_iters, num_obs) {
    Pm <- theta[1]
    Pp <- theta[2]
    Y <- simulate_cell(Yinit, rows, cols, Pm, Pp, sim_iters, num_obs)
    return(Y)
}

#' The function \code{cell_sum(Y,sum_options)} calculates the summary statistics for this example.
#' @rdname cell
#' @export
cell_sum <- function(Y, Yinit) {
    num_obs = dim(Y)[3]
    summ_stat = numeric(num_obs+1)

    # Hamming distances between cell locations across time
    summ_stat[1] = sum(abs(Yinit-Y[, , 1]))
    for (i in 2:num_obs) {
        summ_stat[i] = sum(abs(Y[, , i-1]-Y[, , i]))
    }

    # Total number of cells in the final time period
    summ_stat[num_obs + 1] = sum(Y[, , num_obs])

	return(summ_stat)
}

#' The function \code{cell_prior(theta)} evaluates the prior at the chosen parameters.
#' @rdname cell
#' @export
cell_prior <- function(theta) {
    theta[1] > 0 & theta[1] < 1 & theta[2] > 0 & theta[2] < 1
}
