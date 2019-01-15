#' S3 reference class of the result from tuning to select the optimal penalty for BSLasso
#'
#' @description Two functions (print and plot) are provided for class "penbsl".
#'
#' @param x         A "penbsl" class object, typically the output of function \code{\link{selectPenalty}}.
#' @param logscale  A logical indicator whether the x-axis (penalty) should be log transformed. The
#' default is \code{TRUE}.
#' @param digits    The number of digits to print.
#' @param ...       Other arguments.
#' @inheritParams show.bsl
#'
#' @author 								Ziwen An, Leah F. South and Christopher C. Drovandi
#'
#' @name penbsl
#' @aliases penbsl
NULL

#' print a "penbsl" class result
#' @rdname penbsl
#' @method print penbsl
#' @export
##' @rawNamespace S3method(print,penbsl)
print.penbsl <- function(x, digits = max(3L, getOption("digits") - 4L), ...) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (length(x$resultsDF)) {
        r1 <- x$resultsDF[which(x$resultsDF$sigmaOpt), c('n', 'penalty', 'sigma')]
        cat('Penalty selected based on the standard deviation of the loglikelihood:\n')
        print(format(r1, digits = digits))
    } else {
        cat("No result to show\n")
    }
    return(invisible(r1))
}

#' The function \code{plot.penbsl} can be used to plot the results from tuning to select the optimal penalty for BSLasso.
#' @rdname penbsl
#' @method plot penbsl
#' @export
plot.penbsl <- function(x, logscale = TRUE, ...) {
    sigma <- sigmaOpt <- penalty <- logPenalty <- NULL # to satisfy R CMD check
    x <- x$resultsDF
    a <- floor(sqrt(length(unique(x$n))))
    b <- ceiling(length(unique(x$n)) / a)
    nRepeats <- sapply(unique(x$n), FUN = function(xx) sum(x['n'] == xx))
    yPosSigma <- sapply(unique(x$n), FUN = function(xx) mean(range(x[which(x$n == xx), 'sigma'])))
    textYSigma <- c(unlist(mapply(yPosSigma, nRepeats, FUN = rep)))
    sigmaOpt.NA <- x$sigmaOpt
    sigmaOpt.NA[which(sigmaOpt.NA == FALSE)] <- NA
    sigmaOpt.NA[which(sigmaOpt.NA == TRUE)] <- 1
    logPenalty <- log(x$penalty)
    if (logscale) {
        ggplot(data = x, aes(x = logPenalty, y = sigma)) +
            geom_line(color = 'darkblue', linetype = 'dashed', size = 1) +
            facet_wrap( ~ n, scales = 'free', nrow = a, ncol = b, labeller = label_both) +
            geom_vline(aes(xintercept = logPenalty*sigmaOpt.NA), na.rm = TRUE, color = 'forestgreen', linetype = 4) +
            geom_label(aes(x = logPenalty*sigmaOpt.NA, y = textYSigma, label = paste('lambda == ', round(penalty, 3))),
                       hjust = 0.5, vjust = "inward", parse = TRUE, color = 'white', fill = '#FE66A9', size = 2.7,
                       alpha = 0.8, na.rm = TRUE) +
            labs(x = 'log(penalty)', title = 'Penalty Selection') +
            theme(plot.title = element_text(size = 14, hjust = 0.5)) +
            theme(strip.text.x = element_text(size = 12, face = 'bold'), axis.title = element_text(size = 12))
    } else {
        ggplot(data = x, aes(x = penalty, y = sigma)) +
            geom_line(color = 'darkblue', linetype = 'dashed', size = 1) +
            facet_wrap( ~ n, scales = 'free', nrow = a, ncol = b, labeller = label_both) +
            geom_vline(aes(xintercept = penalty*sigmaOpt.NA), na.rm = TRUE, color = 'forestgreen', linetype = 4) +
            geom_label(aes(x = penalty*sigmaOpt.NA, y = textYSigma, label = paste('lambda == ', round(penalty, 3))),
                       hjust = 0.5, vjust = "inward", parse = TRUE, color = 'white', fill = '#FE66A9', size = 2.7,
                       alpha = 0.8, na.rm = TRUE) +
            labs(x = 'penalty', title = 'Penalty Selection') +
            theme(plot.title = element_text(size = 14, hjust = 0.5)) +
            theme(strip.text.x = element_text(size = 12, face = 'bold'), axis.title = element_text(size = 12))
    }
}
