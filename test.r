analysisofvariance <- function(formula, data) {
    mf <- match.call()
    mf[[1]] <- quote(model.frame)
    data <- eval(mf, envir=parent.frame())
    aov.group.colname <- as.character(attr(attributes(data)$terms, 'term.labels'))
    aov.data.colname <- as.character(formula(attributes(data)$terms)[[2]])
    aov.group.names <- sort(unique(data[[aov.group.colname]]))

    aov.mean <- mean(data[[aov.data.colname]])
    aov.n <- dim(data)[1]
    aov.k <- length(aov.group.names)
    aov.group.mean <- tapply(
        data[[aov.data.colname]],
        data[[aov.group.colname]],
        mean
    )
    aov.group.n <- table(data[[aov.group.colname]])

    aov.tss <- sum((data[[aov.data.colname]] - aov.mean)^2)
    aov.wss <- sum(apply(
        data, 1,
        function(x) {
            (as.numeric(x[aov.data.colname]) - aov.group.mean[x[aov.group.colname]])^2
        }
    ))
    aov.wms <- aov.wss / (aov.n - aov.k)
    aov.bss <- sum(apply(
        data, 1,
        function(x) {
            (as.numeric(aov.group.mean[x[aov.group.colname]]) - aov.mean)^2
        }
    ))
    aov.bms <- aov.bss / (aov.k - 1)
    aov.F <- aov.bms / aov.wms
    aov.F.p <- pf(aov.F,df1=aov.k-1, df2=aov.n-aov.k, lower.tail=FALSE)

    ## calculate group T test
    aov.combine <- combn(seq(aov.k), 2)
    aov.t <- numeric()
    aov.t.p <- numeric()
    for (i in seq(dim(aov.combine)[2])) {
        aname <- aov.group.names[aov.combine[1,i]]
        bname <- aov.group.names[aov.combine[2,i]]
        thename <- paste(aname, bname, sep=':')
        aov.t[thename] <- (
            aov.group.mean[aname] - aov.group.mean[bname]
        ) / sqrt(aov.wms * (1/aov.group.n[aname] + 1/aov.group.n[bname]))
        p1 <- pt(aov.t[thename], df=aov.n-aov.k)
        aov.t.p[thename] <- 2 * min(p1, 1-p1)
    }
    result <- list(
        F=aov.F, F.p=aov.F.p, t=aov.t, t.p=aov.t.p
    )
    return(result)
}
