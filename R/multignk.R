#' The multivariate G&K example
#'
#' @description Here we provide the data and tuning parameters required to reproduce
#' the results from the multivariate G & K (Drovandi and Pettitt, 2011) example from An et al. (2018).
#'
#' @param theta_tilde   A vector with 15 elements for the proposed model parameters.
#' @param T             The number of observations in the data.
#' @param J             The number of variables in the data.
#' @param bound         A matrix of boundaries for the uniform prior.
#' @param y				A \code{T} \eqn{x} \code{J} matrix of data.
#'
#' @details
#' It is not practical to give a reasonable explanation of this example through R documentation
#' given the number of equations involved. We refer the reader to the BSLasso paper (An et al., 2018)
#' at \url{https://doi.org/10.1080/10618600.2018.1537928} for information on the model and summary statistic used in this example.
#'
#' @section An example dataset:
#'
#' We use the foreign currency exchange data available from \url{http://www.rba.gov.au/statistics/historical-data.html}
#' as in An et al. (2018).
#'
#' \itemize{
#'  \item \code{data}:  A \code{1651} \eqn{x} \code{3} matrix of data.
#'  \item \code{sim_options}: Values of \code{sim_options} relevant to this example.
#'  \item \code{start}: A vector of suitable initial values of the parameters for MCMC.
#'  \item \code{cov}: Covariance matrix of the multivariate normal random walk, in the form of a \eqn{15 x 15} matrix.
#' }
#'
#' @examples
#' \donttest{
#' require(doParallel) # You can use a different package to set up the parallel backend
#' require(MASS)
#' require(elliplot)
#'
#' # Loading the data for this example
#' data(mgnk)
#'
#' # Performing BSL (reduce the number of iterations M if desired)
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' resultMgnkBSL <- bsl(mgnk$data, n = 60, M = 80000, theta0 = mgnk$start, covRandWalk = mgnk$cov,
#'     fnSim = mgnk_sim, fnSum = mgnk_sum, simArgs = mgnk$sim_options,
#'     parallel = TRUE, parallelArgs = list(.packages = c('BSL', 'MASS'), .export = 'ninenum'),
#'     thetaNames = expression(a[1],b[1],g[1],k[1],a[2],b[2],g[2],k[2],a[3],b[3],g[3],k[3],
#'         delta[12],delta[13],delta[23]), verbose = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' show(resultMgnkBSL)
#' summary(resultMgnkBSL)
#' plot(resultMgnkBSL, which = 2, thin = 20)
#'
#' # Performing tuning for BSLasso
#' lambda_all <- list(exp(seq(-2.5,0.5,length.out=20)), exp(seq(-2.5,0.5,length.out=20)),
#'                    exp(seq(-4,-0.5,length.out=20)), exp(seq(-5,-2,length.out=20)))
#'
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' sp_mgnk <- selectPenalty(ssy = mgnk_sum(mgnk$data), n = c(15, 20, 30, 50), lambda_all,
#'     theta = mgnk$start, M = 100, sigma = 1.5, fnSim = mgnk_sim,
#'     fnSum = mgnk_sum, simArgs = mgnk$sim_options, standardise = TRUE,
#'     parallelSim = TRUE, parallelSimArgs = list(.packages = c('BSL', 'MASS'), .export = 'ninenum'),
#'     parallelMain = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' sp_mgnk
#' plot(sp_mgnk)
#'
#' # Performing BSLasso with a fixed penalty (reduce the number of iterations M if desired)
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' resultMgnkBSLasso <- bsl(mgnk$data, n = 20, M = 80000, theta0 = mgnk$start, covRandWalk = mgnk$cov,
#'      fnSim = mgnk_sim, fnSum = mgnk_sum, simArgs = mgnk$sim_options,
#'      penalty = 0.3, standardise = TRUE, parallel = TRUE,
#'      parallelArgs = list(.packages = c('BSL', 'MASS'), .export = 'ninenum'),
#'      thetaNames = expression(a[1],b[1],g[1],k[1],a[2],b[2],g[2],k[2],a[3],b[3],g[3],k[3],
#'          delta[12],delta[13],delta[23]), verbose = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' show(resultMgnkBSLasso)
#' summary(resultMgnkBSLasso)
#' plot(resultMgnkBSLasso, which = 2, thin = 20)
#'
#'
#' # Performing semiBSL (reduce the number of iterations M if desired)
#' # Opening up the parallel pools using doParallel
#' cl <- makeCluster(detectCores() - 1)
#' registerDoParallel(cl)
#' resultMgnkSemiBSL <- bsl(mgnk$data, n = 60, M = 80000, theta0 = mgnk$start, covRandWalk = mgnk$cov,
#'     fnSim = mgnk_sim, fnSum = mgnk_sum, simArgs = mgnk$sim_options, method = 'semiBSL',
#'     parallel = TRUE, parallelArgs = list(.packages = c('BSL', 'MASS'), .export = 'ninenum'),
#'     thetaNames = expression(a[1],b[1],g[1],k[1],a[2],b[2],g[2],k[2],a[3],b[3],g[3],k[3],
#'         delta[12],delta[13],delta[23]), verbose = TRUE)
#' stopCluster(cl)
#' registerDoSEQ()
#' show(resultMgnkSemiBSL)
#' summary(resultMgnkSemiBSL)
#' plot(resultMgnkSemiBSL, which = 2, thin = 20)
#'
#' # Plotting the results together for comparison
#' # plot using the R default plot function
#' par(mar = c(4, 4, 1, 1), oma = c(0, 1, 2, 0))
#' combinePlotsBSL(list(resultMgnkBSL, resultMgnkBSLasso, resultMgnkSemiBSL), which = 1, thin = 20,
#'     label = c('bsl', 'bslasso', 'semiBSL'), col = c('red', 'blue', 'green'), lty = 2:4, lwd = 1)
#' mtext('Approximate Univariate Posteriors', outer = TRUE, line = 0.75, cex = 1.2)
#'
#' # plot using the ggplot2 package
#' library("grid") # use function "unit" from package "grid" to fix margins, remove if wanted
#' combinePlotsBSL(list(resultMgnkBSL, resultMgnkBSLasso, resultMgnkSemiBSL), which = 2, thin = 20,
#'     label=c('bsl   ','bslasso   ','semiBSL'), options.color=list(values=c('red','blue','green')),
#'     options.linetype = list(values = 2:4), options.size = list(values = rep(1, 3)),
#'     options.theme = list(plot.margin = unit(rep(0.03,4),"npc"), axis.title = element_text(size=12),
#'         axis.text = element_text(size = 8), legend.text = element_text(size = 12)))
#' }
#'
#' @references
#' An, Z., South, L. F., Nott, D. J. &  Drovandi, C. C. (2018). Accelerating Bayesian synthetic
#' likelihood with the graphical lasso. Journal of Computational and Graphical Statistics.
#' \url{https://doi.org/10.1080/10618600.2018.1537928}
#'
#' Drovandi, C. C. and Pettitt, A. N. (2011). Likelihood-free Bayesian estimation of multivariate
#' quantile distributions. Computational Statistics and Data Analysis, 55(9):2541-2556.
#'
#' @author 								Ziwen An, Leah F. South and Christopher C. Drovandi
#'
#' @name mgnk
#' @docType data
NULL
## NULL

# quantile function of a g-and-k distribution
qgnk <- function(z, a, b, g, k) {
    e <- exp(- g * z)
    a + b * (1 + 0.8 * (1 - e) / (1  + e)) * (1 + z^2) ^ k * z
}

logTransform <- function(x, bound) {
    x_tilde <- numeric(4)
    for (i in 1 : 4) {
        x_tilde[i] <- log((x[i] - bound[i, 1]) / (bound[i, 2] - x[i]))
    }
    return(x_tilde)
}

backLogTransform <- function(x_tilde, bound) {
    x <- numeric(4)
    for (i in 1 : 4) {
        x[i] <- (bound[i, 1] + bound[i, 2] * exp(x_tilde[i])) / (1 + exp(x_tilde[i]))
    }
    return(x)
}

reparaCorr <- function(theta_corr, J) {
    Sigma <- diag(J)
    count <- 1

    for (i in 1 : (J-1)) {
        for (j in (i+1) : J) {
            Sigma[i, j] <- Sigma[j, i] <- theta_corr[count]
            count <- count + 1
        }
    }

    L <- t(chol(Sigma))
    gamma <- matrix(0, J, J)
    w <- numeric(choose(J, 2))
    count <- 1

    for (i in 2 : J) {
        gamma[i, 1] <- acos(L[i, 1])
    }
    for (j in 2 : (J-1)) {
        for (i in (j+1): J) {
            gamma[i, j] <- acos((L[i, j]) / (prod(sin(gamma[i, 1:(j-1)]))))
        }
    }

    for (i in 2 : J) {
        for (j in 1: (i-1)) {
            w[count] <- log(gamma[i, j] / (pi - gamma[i, j]))
            count <- count + 1
        }
    }

    return(list(w = w, Sigma = Sigma))
}

backReparaCorr <- function(w, J) {
    G <- array(0, c(J, J))
    count <- 1

    for (i in 2 : J) {
        for (j in 1 : (i-1)) {
            G[i, j] <- pi / (1 + exp(-w[count]))
            count <- count + 1
        }
    }

    L <- array(0, c(J, J))
    L[1, 1] <- 1
    for (i in 2 : J) {
        L[i, 1] <- cos(G[i, 1])
        L[i, i] <- prod(sin(G[i, 1 : (i-1)]))
    }
    for (i in 3 : J) {
        for (j in 2 : (i-1)) {
            L[i, j] <- prod(sin(G[i, 1 : (j-1)])) * cos(G[i, j])
        }
    }

    Sigma <- L %*% t(L)
    theta_corr <- numeric(choose(J, 2))
    count <- 1

    for (i in 1 : (J - 1)) {
        for (j in (i + 1) : J) {
            theta_corr[count] <- Sigma[i, j]
            count <- count + 1
        }
    }

    return(theta_corr)

}

paraTransformGnk <- function(theta, J, bound) {
    if (J == 1L) {
        theta_tilde <- logTransform(theta, bound)
    } else {
        theta_tilde <- numeric(length(theta))
        for (i in 1:J) {
            theta_tilde[(4*i-3) : (4*i)] <- logTransform(theta[(4*i-3) : (4*i)], bound)
        }
        theta_tilde[(4*J + 1) : length(theta_tilde)] <- reparaCorr(tail(theta, -4*J), J)$w
    }
    return(theta_tilde)
}

paraBackTransformGnk <- function(theta_tilde, J, bound) {
    if (J == 1L) {
        theta <- backLogTransform(theta_tilde, bound)
    } else {
        theta <- numeric(length(theta_tilde))
        for (i in 1:J) {
            theta[(4*i-3) : (4*i)] <- backLogTransform(theta_tilde[(4*i-3) : (4*i)], bound)
        }
        theta[(4*J + 1) : length(theta)] <- backReparaCorr(tail(theta_tilde, -4*J), J)
    }
    return(theta)
}

#' The function \code{mgnk_sim} simulates from the multivariate G & K model.
#' @rdname mgnk
#' @export
mgnk_sim <- function(theta_tilde, T, J, bound) {
    theta <- paraBackTransformGnk(theta_tilde, J, bound)
    if (J == 1) {
        theta_gnk <- theta
        Sigma <- 1
    } else {
        theta_gnk <- head(theta, 4*J)
        theta_corr <- tail(theta, -4*J)

        if (length(theta_corr) != choose(J, 2)) {
            stop('wrong parameter length or dimension')
        }

        Sigma <- reparaCorr(theta_corr, J)$Sigma
    }

    y <- array(0, c(T, J))
    zu <- mvrnorm(n = T, mu = numeric(J), Sigma = Sigma)

    for (i in 1 : J) {
        y[, i] <- qgnk(zu[, i], theta[4*(i-1) + 1], theta[4*(i-1) + 2], theta[4*(i-1) + 3], theta[4*(i-1) + 4])
    }
    return(y)
}

summStatRobust <- function(x) {
    T <- length(x)
    ssx <- numeric(4)
    octile <- elliplot::ninenum(x)[2:8]
    ssx[1] <- octile[2]
    ssx[2] = octile[3] - octile[1]
    ssx[3] = (octile[3] + octile[1] - 2*octile[2]) / ssx[2]
    ssx[4] = (octile[7] - octile[5] + octile[3] - octile[1]) / ssx[2]
    return(ssx)
}

normScore <- function(x, y) {
    n <- length(x)
    r0 <- 1 : n
    z1 <- qnorm(rank(x) / (n + 1))
    z2 <- qnorm(rank(y) / (n + 1))
    c <- qnorm(r0 / (n + 1))
    r <- sum(z1 * z2) / sum(c ^ 2)
    norm_score <- 0.5 * log((1 + r) / (1 - r))
    return(norm_score)
}

#' The function \code{mgnk_sum(y)} calculates the summary statistics for the multivariate G & K example.
#' @rdname mgnk
#' @export
mgnk_sum <- function(y) {
    J <- ncol(y)
    ssxRobust <- c(apply(y, MARGIN = 2, FUN = summStatRobust))
    if (J == 1L) {
        return(ssx = ssxRobust)
    } else {
        ssxNormScore <- numeric(choose(J, 2))
        count <- 1
        for (i in 1 : (J - 1)) {
            for (j in (i + 1) : J) {
                ssxNormScore[count] <- normScore(y[, i], y[, j])
                count <- count + 1
            }
        }
        return(ssx = c(ssxRobust, ssxNormScore))
    }
}
