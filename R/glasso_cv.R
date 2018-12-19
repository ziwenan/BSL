glasso_cv <- function(ssx, penalty, K = 3L) {
    n <- nrow(ssx)
    p <- ncol(ssx)
    cv_unit <- numeric(K)
    folds <- cvFolds(n, K = K)
    for (k in 1 : K) {
        nk <- length(folds$which == k)
        ssx_curr <- ssx[folds$which != k, ]
        S_curr <- cov(ssx[folds$which == k, ])
        gl <- glasso(cov(ssx_curr), penalty)
        Omega <- gl$wi
        cv_unit[k] <- -nk * log(det(Omega)) + nk * sum(diag(Omega %*% S_curr))
    }
    cv <- sum(cv_unit)
    return(cv)
}
