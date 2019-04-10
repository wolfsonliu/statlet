# -*- coding: utf-8 -*-

#' Create rank matrix
#'
#' Convert a set of ranked lists into a rank matrix
#'
#' The lists are converted to a format that is used by aggregateRanks. If partial
#' rankings are given to the function, all the missing values are subtituted by the
#' maximum rank N, which can be specified manually. This parameter has a very strong
#' influence on the performance of RRA algorithm, therfore it should be reasonably
#' accurate. If the N is different for the gene lists, it can be also given as a vector.
#'
#' Parameter full is used, when full rankings are given, but the sets of ranked elements
#' do not match perfectly. Then the structurally missing values are substituted with
#' NA-s.
#'
#' @param x list of preference lists
#' @return A matrix, with as many columns as input rankings and rows as unique elements
#' in all the rankings combined.
#' @author  Raivo Kolde \email{rkolde@@gmail.com}
#' #' @examples
#' # Make sample input data
#' glist <- list(sample(letters, 4), sample(letters, 10), sample(letters, 12))
#'
#' r = rankMatrix(glist)
#' r = rankMatrix(glist, full = TRUE)
#'
#' # Use real data
#' data(cellCycleKO)
#' r = rankMatrix(cellCycleKO$gl, N = cellCycleKO$N)
#'
#' @export

rankMatrix <- function(x, ...) {
    UseMethod('rankMatrix', x)
}


#' @inherit rankMatrix
#' @param N number of all rankable elements
#' @param full logical showing if the given rankings are complete
#' @author  Wolfson Liu \email{wolfsonliu@@gmail.com}
#' @export

rankMatrix.default <-  function(x, N=NA, full=FALSE, ...) {
    u <- unique(c(x, recursive=TRUE))
    if (all(is.na(N))) {
        N <- length(u)
    }
    if (!full) {
        rmat <- matrix(
            1, nrow=length(u), ncol=length(x),
            dimnames=list(u, names(x))
        )
        if (length(N) == 1) {
            N <- rep(N, ncol(rmat))
        }
    } else {
        rmat <- matrix(
            NA, nrow=length(u), ncol=length(x),
            dimnames=list(u, names(x))
        )
        N <- unlist(lapply(x, length))
    }
    for (i in 1:length(x)) {
        rmat[x[[i]], i] <- (1:length(x[[i]])) / N[i]
    }
    return(rmat)
}

#' @inherit rankMatrix
#' @param direction c('low', 'high') output the rmat with gene of decrease or increase in experiment
#' @param group c('gene', 'sgrna') consider gene level or sgrna level
#' @export

rankMatrix.matrix <- function(x,
                              normalization=FALSE,
                              byrow=TRUE,
                              na.last=TRUE,
                              ties.method=c(
                                  'average', 'first', 'last',
                                  'random', 'max', 'min'
                              ),
                              ...) {
    ties.method <- match.arg(ties.method)

    if (byrow) {
        x <- t(x)
    } else {}

    if (is.null(colnames(x))) {
        cols <- seq(dim(x)[2])
    } else {
        cols <- colnames(x)
    }

    if (is.null(rownames(x))) {
        rows <- seq(dim(x)[1])
    } else {
        rows <- rownames(x)
    }

    rk <- matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])

    for (i in cols) {
        rk[,i] <- rank(
            x[, i],
            na.last=na.last,
            ties.method=ties.method
        )
    }

    if (normalization) {
        rk <- rk / length(rows)
    } else {}

    if (byrow) {
        rk <- t(rk)
    } else {}

    return(rk)
}

rankMatrix.data.frame <- rankMatrix.matrix
## -----------------
