#' @include S4Functions.R
NULL

setOldClass('difftime')

setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("functionOrNULL", c("function", "NULL"))

#' An S4 class to represent a bsl result.
#' @export
setClass("bsl", slots = c(theta = "matrix", loglike = "numeric", loglikeAll = "numeric",
    acceptanceRate = "numeric", earlyRejectionRate = "numeric", call = "call",
    y = "ANY", n = "numeric", M = "numeric", theta0 = "numeric", covRandWalk = "matrix",
    fnSim = "function", fnSum = "function", method = "character", shrinkage = "characterOrNULL",
    penalty = "numericOrNULL", fnPrior = "functionOrNULL", simArgs = "listOrNULL", sumArgs = "listOrNULL",
    logitTransform = "logical", logitTransformBound = "matrixOrNULL", standardise = "logical",
    parallel = "logical", parallelArgs = "listOrNULL", thetaNames = "expression",
    time = "difftime"))

setMethod("initialize", "bsl", initialize.bsl)
setValidity("bsl", check.bsl)


#' @export
#setGeneric('show')
setMethod('show', signature(object = 'bsl'), show.bsl)
#setGeneric('show', show.bsl)

#' @export
setGeneric('summary')
setMethod('summary', 'bsl', summary.bsl)

#' @export
setGeneric('plot')
setMethod('plot', signature = c(x = 'bsl', y = 'missing'), plot.bsl)
