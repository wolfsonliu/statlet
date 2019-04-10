# -*- coding: utf-8 -*-

rowVars <- function(x, na.rm=FALSE, ...)
{
    n <- rowSums(!is.na(x))
    n[n <= 1] <- NA
    return(
        rowSums(
            (x - rowMeans(x,na.rm=na.rm, ...))^2,
            na.rm=na.rm
        )/(n - 1)
    )
}

rowMedians <- function(x, na.rm=FALSE, ...)
{
    apply(x, MARGIN=1, median, na.rm=na.rm, ...)
}

colVars <- function(x, na.rm=FALSE, ...)
{
    n <- colSums(!is.na(x))
    n[n <= 1] <- NA
    return(
        rowSums(
            (t(x) - colMeans(x, na.rm=na.rm, ...))^2,
            na.rm=na.rm
        )/(n-1)
    )
}

colMedians <- function(x, na.rm=FALSE, ...)
{
    apply(x, MARGIN=2, median, na.rm=na.rm, ...)
}


geommean <- function(x, na.rm=FALSE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / sum(!is.na(x)))
}

#' @title Normalize data
#'
#' The total sequence number can affect the sgRNA counts.
#' The data should be normalized before investigation.
#'
#' @param x a vector or matrix, containing numeric, integer.
#' @param method method used in normalization.
#'     rpm, z, mean, median could be used for vector or matrix.
#'     UQS, quantile could be used for matrix.
#'
#' @examples
#' x <- seq(1:10)
#'
#' y <- cbind(seq(1,10), seq(2,20,2), 3)
#'
#' normalize(x, method='rpm')
#'
#' normalize(y, 'rpm', na.rm=TRUE, col.as.rep=FALSE)
#'
#' @export
normalize <- function(x, method, ...)
{
    UseMethod('normalize', x)
}

#' @param sizefactor number, used in rpm normalization
#'
#' @details
#' rpm: reads per million, sizefactor usually set to 1 million
#' z: z score
#' mean: all counts divided by mean
#' median: all counts divided by median
#'
#' @describeIn normalize
#' @export
normalize.vector <- function(x,
                             method=c(
                                 'rpm', 'z', 'mean', 'median'
                             ),
                             ...,
                             na.rm=FALSE,
                             sizefactor=10^6) {
    method <- match.arg(method)
    if (length(x) == 1) {
        stop('The lenght of x should be larger than 1.')
    } else {}

    if (na.rm) {
        x <- x[!is.na(x)]
    } else {}

    if (method == 'rpm') {
        normx <- x / sum(x, na.rm=na.rm) * sizefactor
    } else if (method == 'z') {
        normx <- (x - mean(x, na.rm=na.rm)) / sqrt(var(x, na.rm=na.rm))
    } else if (method == 'mean') {
        normx <- x / mean(x, na.rm=na.rm)
    } else if (method == 'median') {
        normx <- x / median(x, na.rm=na.rm)
    }

    normx
}

#' @details
#' UQS: upper quantile scaling
#' quantile: quantile normalization
#'
#' @describeIn normalize
#' @export
normalize.matrix <- function(x,
                             method=c(
                                 'rpm', 'z', 'mean', 'median',
                                 'UQS', 'quantile'
                             ),
                             ...,
                             row.as.rep=FALSE,
                             na.fill=NA,
                             na.rm=TRUE,
                             sizefactor=10^6) {
    method <- match.arg(method)

    if (!row.as.rep) {
        rx <- t(x)
    } else {
        rx <- x
    }
    ## x dim to replicant, item

    rx[is.na(rx)] <- na.fill

    if (method == 'rpm') {
        sums <- rowSums(rx, na.rm=na.rm) # replicants num length
        normx <- rx / sums * sizefactor  # broadcast sums
    } else if (method == 'z') {
        means <- rowMeans(rx, na.rm=na.rm)
        vars <- rowVars(rx, na.rm=na.rm)

        if (any(vars)) {
            stop('variance is zero')
        } else {}

        normx <- (rx - means) / sqrt(vars)
    } else if (method == 'mean') {
        means <- rowMeans(rx, na.rm=na.rm)

        normx <- rx / means              # broadcast means
    } else if (method == 'median') {
        medians <- rowMedians(rx, na.rm=na.rm)

        normx <- rx / medians            # broadcast medians
    } else if (method == 'UQS') {
        ## upper quantile scale
        uq <- apply(rx, MARGIN=1, fivenum)[4, ] # replicants num length
        uqgm <- geomean(uq)                    # geometric mean of uq
        normx <- rx / uq * uqgm          # broadcast upper quantile
    } else if (method == 'quantile') {
        rx.sorted <- apply(rx, MARGIN=1, sort, decreasing=FALSE)
        repmean <- sort(rowMeans(rx.sorted, na.rm=TRUE))
        rx.rank <- t(apply(rx,1,rank,ties.method="min"))
        normx <- matrix(repmean[rx.rank], nrow=dim(rx.rank)[1])
    }

    if (!row.as.rep) {
        normx <- t(normx)
    } else {
    }
    colnames(normx) <- colnames(x)
    rownames(normx) <- rownames(x)
    return(normx)
}

#' @details
#' UQS: upper quantile scaling
#' quantile: quantile normalization
#'
#' @describeIn normalize
#' @export
normalize.data.frame <- normalize.matrix
## -----------------
