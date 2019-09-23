semiparaKernelEstimate2 <- function (ssy, ssx, kernel = 'gaussian', shrinkage = NULL,
                                    penalty = NULL, log = TRUE) {
    if (!is.null(shrinkage)) {
        flagShrinkage <- TRUE
        shrinkage <- match.arg(shrinkage, c("glasso", "Warton"))
    } else {
        flagShrinkage <- FALSE
    }
    if (!flagShrinkage && !is.null(penalty)) {
        warning('"penalty" will be ignored since no shrinkage method is specified')
    }
    if (flagShrinkage && is.null(penalty)) {
        stop('a penalty value must be specified to enable shrinkage estimation')
    }

    n <- nrow(ssx)
    ns <- ncol(ssx)

    pdfy <- yu <- numeric(ns)
    for (j in 1 : ns) {
        f <- density(ssx[, j], kernel = kernel, n = 512, from = min(ssx[,j],ssy[j]),
                     to = max(ssx[,j],ssy[j]))
        approxy <- approx(f$x, f$y, ssy[j])
        pdfy[j] <- approxy$y
        yu[j] <- mean(kernelCDF((ssy[j] - ssx[, j]) / f$bw, kernel))
    }

    # Gaussian rank correlation
    rhohat <- gaussianRankCorr(ssx, TRUE)

    if (!is.null(shrinkage)) {
        RHOHAT <- p2P(rhohat)
        Sigma <- switch(shrinkage,
                        'glasso' = glasso(RHOHAT, rho = penalty, penalize.diagonal = FALSE)$w,
                        'Warton' = corrWarton(RHOHAT, penalty))
        rhohat <- P2p(Sigma)
    }

    # density
    copula <- normalCopula(rhohat, dim = ns, dispstr = 'un')
    if (log) {
        f <- dCopula(yu, copula, log = TRUE) + sum(log(pdfy))
    } else {
        f <- dCopula(yu, copula, log = FALSE) * prod(pdfy)
    }

    return(f)
}
