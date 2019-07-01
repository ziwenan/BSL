#' Gaussian rank correlation
#' @description This function computes the Gaussian rank correlation of Boudt et al. (2012).
#' @param x         A numeric matrix representing data where the number of rows is the number of 
#' independent data points and the number of columns are the number of variables in the dataset.
#' @param vec      A logical argument indicating if the vector of correlations should be returned 
#' instead of a matrix.
#' @return         Gaussian rank correlation matrix (default) or a vector of pair correlations.
#' @references
#' Boudt, K., Cornelissen, J., and Croux, C. (2012). The Gaussian rank correlation estimator: 
#' robustness properties. Statistics and Computing, 22(2):471-483.
#' @rdname gaussianRankCorr
#' @examples
#' data(ma2)
#' set.seed(100)
#' 
#' # generate 1000 simualtions from the ma2 simulation function
#' x <- t(replicate(1000, ma2_sim(ma2$start, 10)))
#' 
#' corr1 <- cor(x) # traditional correlation matrix
#' corr2 <- gaussianRankCorr(x) # Gaussian rank correlation matrix
#' par(mfrow = c(1, 2))
#' image(corr1, main = 'traditional correlation matrix')
#' image(corr2, main = 'Gaussian rank correlation matrix')
#' 
#' std <- apply(x, MARGIN = 2, FUN = sd) # standard deviations
#' cor2cov(gaussianRankCorr(x), std) # convert to covariance matrix
#' 
#' @seealso       \code{\link{cor2cov}} for converting a correlation matrix to a covariance matrix.
#' @export
gaussianRankCorr <- function(x, vec = FALSE) {
    n <- nrow(x)
    p <- ncol(x)
    r <- apply(x, FUN = rank, MARGIN = 2, ties.method = 'average')
    rqnorm <- qnorm(r / (n + 1))
    den <- sum((qnorm(1 : n / (n + 1))) ^ 2)
    res <- unlist(sapply(1:(p-1), FUN = function(i) c(rqnorm[, i] %*% rqnorm[, (i+1):p] / den)))
    if (!vec) {
	   res <- p2P(res)
	}
	return (res)
}

#' Convert a correlation matrix to a covariance matrix
#' @description This function converts a correlation matrix to a covariance matrix
#' @param corr    The correlation matrix to be converted. This must be symmetric.
#' @param std     A vector that contains the standard deviations of the variables in the correlation matrix.
#' @return        The covariance matrix.
#' @export
cor2cov <- function(corr, std) {
	outer(std, std) * corr
}