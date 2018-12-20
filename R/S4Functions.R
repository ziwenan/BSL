initialize.bsl <- function(.Object, theta = numeric(0), M = numeric(0), time = numeric(0), ...) {
    .Object <- callNextMethod()
    if (length(theta) == 0) {
        .Object@theta <- matrix(0, 10, 1)
    }
    if (length(M) == 0) {
        .Object@M <- 10
    }
    if (length(time) == 0) {
        .Object@time <- as.difftime(0, units = 'secs')
    }
    validObject(.Object)
    return (.Object)
}



check.bsl <- function(object) {
    if (any(length(object@theta) == 0, length(object@M) == 0)) { # slots that must include in bsl class
        warnings('empty slot "theta" or "M" in the "bsl" object')
    } else {
        errors <- character()
        p <- ncol(object@theta)
        M <- nrow(object@theta)
        if (M != object@M) {
            msg <- paste('The number of rows of theta', M, 'does not match the number of iterations M', object@M)
            error <- c(errors, msg)
        }

        temp <- length(object@loglike)
        if (temp != 0 && temp != M) {
            msg <- paste('The number of iterations M', M, 'does not match the length of loglike', temp)
            error <- c(errors, msg)
        }

        temp <- length(object@loglikeAll)
        if (temp != 0 && temp != M) {
            msg <- paste('The number of iterations M', M, 'does not match the length of loglikeAll', temp)
            error <- c(errors, msg)
        }

        temp <- length(object@theta0)
        if (temp != 0 && temp != p) {
            msg <- paste('The length of theta0', temp, 'does not match the number of parameters', p)
            error <- c(errors, msg)
        }

        if (nrow(object@covRandWalk) != p || ncol(object@covRandWalk) != p) {
            msg <- paste('covRandWalk must be a', p, 'by', p, 'square matrix')
            error <- c(errors, msg)
        }

        if (!is.null(object@logitTransformBound)) {
            if (nrow(object@logitTransformBound) != p || ncol(object@logitTransformBound) != 2L) {
                msg <- paste('logitTransformBound must be a', p, 'by', 2, 'matrix')
                error <- c(errors, msg)
            }
        }

        temp <- length(object@thetaNames)
        if (temp != 0 && temp != p) {
            msg <- paste('The length of thetaNames', temp, 'does not match the number of parameters', p)
            error <- c(errors, msg)
        }

        if (length(errors) == 0) {
            return (TRUE)
        } else {
            return (errors)
        }
    }
}



#' Show method for class "bsl". Display the basic information of a bsl object.
#' @description Display the basic information of a bsl object.
#' @param object   A "bsl" class object to be displayed.
show.bsl <- function(object) {
    digits = max(3L, getOption("digits") - 3L)
    cat("\nCall:\n", paste(deparse(object@call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (nrow(object@theta)) {
        cat("Summary of theta:\n")
        summ <- summary(object@theta)
        attr(summ, 'dimnames') = list(NULL, object@thetaNames)
        print.default(format(summ, digits = digits), print.gap = 2L,
                      quote = FALSE)
    }
    else cat("No theta\n")
    if (length(object@loglike)) {
        cat("Summary of loglikelihood:\n")
        summ <- summary(object@loglike)
        print.default(format(summ, digits = digits), print.gap = 2L,
                      quote = FALSE)
    }
    else cat("No loglikelihood\n")
    if (length(object@acceptanceRate)) {
        cat("Acceptance Rate:\n")
        print.default(format(object@acceptanceRate, digits = digits), print.gap = 2L,
                      quote = FALSE)
    }
    else cat("No acceptance rate\n")
    if (length(object@earlyRejectionRate)) {
        cat("Early Rejection Rate:\n")
        print.default(format(object@earlyRejectionRate, digits = digits), print.gap = 2L,
                      quote = FALSE)
    }
    else cat("No early rejection rate\n")
    cat("\n")
}


#' Summary method for class "bsl"
#' @description Summarise a bsl class object. The effective sample size is scaled with the total number of simulations.
#' @param object     A "bsl" class object to be summarised.
#' @param thetaNames Parameter names to be shown in the summary table. Parameter names of the bsl object will be used by default.
#' @param scale      A constant value to scale up the effective sample size. Default is 1 million.
#' @return A vector of the number of simulations per iteration, acceptance rate of the Markov chain annd scaled effective sample size for each parameter.
summary.bsl <- function(object, thetaNames = NULL, scale = 1000000) {
    theta <- object@theta
    n <- object@n
    M <- nrow(theta)
    p <- ncol(theta)
    if (is.null(thetaNames)) {
        if (!is.null(object@thetaNames)) {
            thetaNames <- object@thetaNames
        } else {
            thetaNames <- vector('expression',p)
            for (i in 1:p) {
                thetaNames[i] <- as.expression(substitute(theta[j],list(j=i)))
            }
        }
    }
    if (length(thetaNames) != p) {
        warning('length of "thetaNames" does not match number of parameters\n')
    }
    if (!is.null(object@acceptanceRate)) {
        accRate = round(object@acceptanceRate, 2)
    }
    else {
        accRate = mean(diff(theta[, 1]) !=0)
    }
    ess <- round(effectiveSize(theta) / n / M * scale, 0)
    summ <- c(n, accRate*100, ess)
    names(summ) <- c('n', 'acc. rate (\\%)', paste('EES', thetaNames))
    return(summ)
}



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
#'
#' @examples
#' \dontrun{
#' # pretend we had a bsl result
#' result <- new('bsl')
#' result@theta <- MASS::mvrnorm(10000, c(0.6, 0.2), diag(c(1, 1)))
#' result@M <- 10000
#'
#' # plot using the R default plot function
#' par(mar = c(5, 4, 1, 2), oma = c(0, 1, 3, 0))
#' plot(result, which = 1, thin = 10, thetaTrue = c(0.6, 0.2),
#'      options.plot = list(cex.main = 1, col = 'red', lty = 2, lwd = 2, main = NA))
#' mtext('Approximate Univariate Posteriors', outer = TRUE, cex = 1.5)
#'
#' # plot using the ggplot2 package
#' plot(result, which = 2, thin = 10, thetaTrue = c(0.6, 0.2),
#'      options.density = list(colour = 'darkblue', fill = 'grey80', size = 1),
#'      options.theme = list(plot.margin = unit(rep(0.05,4), "npc"),
#'                           axis.text = ggplot2::element_text(size = 10)))
#' }
#'
#' @seealso 							\code{\link[BSL:combinePlotsBSL]{combinePlotBSL}} for a function to plot multiple BSL densities.
#' @author 								Ziwen An, Christopher C. Drovandi and Leah F. South
#' @name plot
#' @aliases plot
plot.bsl <- function(x, which = 1L, thin = 1, thetaTrue = NULL, options.plot = NULL,
                     top = 'Approximate Univariate Posteriors', options.density = list(), options.theme = list()) {
    if (which == 1L) {
        if (length(options.density) != 0 || length(options.theme) != 0) {
            warning('"options.density" and "options.theme" are ignored when which = 1')
        }
        marginalPostDefault(x, thin, thetaTrue, options.plot)
    } else if (which == 2L) {
        if (!is.null(options.plot)) {
            warning('"options.plot" is ignored when which = 2')
        }
        marginalPostGgplot(x, thin, thetaTrue, top, options.density, options.theme)
    } else {
        stop('Indicate a supported plot number, 1 for R default density plot or 2 for ggplot density plot')
    }
}


#' Plot the univariate marginal posterior plot of a bsl class object using the R default plot function.
#' @rdname plot
marginalPostDefault <- function(x, thin = 1, thetaTrue = NULL, options.plot = NULL) {
    n <- nrow(x@theta)
    p <- ncol(x@theta)
    a <- floor(sqrt(p))
    b <- ceiling(p / a)
    if (!is.null(thetaTrue) & length(thetaTrue) != p) {
        stop('Length of thetaTrue does not match the number of parameters.')
    }
    thetaNames <- x@thetaNames
    par(mfrow = c(a, b))
    for(k in 1:p) {
        idx <- seq(1, n, thin)
        d <- density(x@theta[idx, k])
        if ('main' %in% names(options.plot)) {
            do.call(plot, c(list(d, xlab = thetaNames[k]), options.plot))
        } else {
            do.call(plot, c(list(d, xlab = thetaNames[k], main = NA), options.plot))
        }
        if (!is.null(thetaTrue)) {
            abline(v = thetaTrue[k], col = 'forestgreen', lty = 3)
        }
    }
    par(mfrow = c(1,1))
}


#' Plot the univariate marginal posterior plot of a bsl class object using the ggplot2 package.
#' @rdname plot
marginalPostGgplot <- function(x, thin = 1, thetaTrue = NULL, top = 'Approximate Univariate Posteriors', options.density = list(), options.theme = list()) {
    n <- nrow(x@theta)
    p <- ncol(x@theta)
    a <- floor(sqrt(p))
    b <- ceiling(p / a)
    if (!is.null(thetaTrue) & length(thetaTrue) != p) {
        stop('Length of thetaTrue does not match the number of parameters.')
    }
    samples <- data.frame(x@theta[seq(1, n, by = thin), ])
    thetaNames <- x@thetaNames
    plist <- list()
    for (i in 1 : p) {
        plist[[i]] <- ggplot(samples, aes_string(x = colnames(samples)[i])) +
            do.call(geom_density, options.density) +
            geom_hline(yintercept = 0, colour = "grey", size = 0.75) + {
                if (!is.null(thetaTrue)) {
                    geom_vline(xintercept = thetaTrue[i], color = 'forestgreen', linetype = 'dashed', size = 0.5)
                }
            } +
            labs(x = thetaNames[i], y = 'density') +
            do.call(theme, options.theme)
    }
    do.call('grid.arrange', c(plist, nrow = a, ncol = b, top = top))
}
