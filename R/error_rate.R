error.rate <- function (clust1, clust2)
{
    clust1 <- unclass(as.ordered(clust1))
    clust2 <- unclass(as.ordered(clust2))
    if ((n = length(clust1)) != length(clust2)) {
        warning("error: length not equal")
        return
    }
    if ((g = length(table(clust1))) != length(table(clust2))) {
        warning("the number of clusters are not equal")
        return
    }
    permute <- function(a) {
        n <- length(a)
        if (n == 1)
            f <- a
        else {
            nm <- gamma(n)
            f <- array(0, c(n, n * nm))
            j <- 1
            for (i in a) {
                f[1, (j - 1) * nm + 1:nm] <- i
                f[-1, (j - 1) * nm + 1:nm] <- permute(setdiff(a,
                  i))
                j <- j + 1
            }
        }
        f
    }
    id <- 1:n
    cmb <- permute(1:g)
    nperm <- ncol(cmb)
    rate <- rep(0, nperm)
    for (i in 1:nperm) {
        tmp <- rep(0, g)
        tc <- rep(0, n)
        for (j in 1:g) tc[clust2 == j] = cmb[j, i]
        for (j in 1:g) {
            tmp1 <- 0
            for (k in (1:g)[-j]) tmp1 <- tmp1 + length(intersect(id[clust1 ==
                j], id[tc == k]))
            tmp[j] <- tmp1
        }
        rate[i] <- sum(tmp)/n
    }
    min(rate)
}


rand.index <- function (LabelA, LabelB)
{
    u <- unclass(as.ordered(LabelA))
    v <- unclass(as.ordered(LabelB))
    if ((N <- length(u)) != length(v))
        stop("Label A and B does not match!")
    row <- max(u)
    col <- max(v)
    nvect <- array(0, c(row, col))
    for (i in 1:row) {
        for (j in 1:col) {
            nvect[i, j] <- sum(u == i & v == j)
        }
    }
    SumsA <- rowSums(nvect)
    SumsB <- colSums(nvect)
    a = 0
    for (i in 1:row) a = a + choose(SumsA[i], 2)
    b = 0
    for (j in 1:col) b = b + choose(SumsB[j], 2)
    c <- a * b/choose(N, 2)
    d = 0
    for (i in 1:row) {
        for (j in 1:col) {
            d = d + choose(nvect[i, j], 2)
        }
    }
    arj <- (d - c)/((a + b)/2 - c)
    a = d
    b = 0
    for (l in 1:row) {
        for (i in 1:(col - 1)) {
            for (j in (i + 1):col) {
                b = b + nvect[l, i] * nvect[l, j]
            }
        }
    }
    c = 0
    for (l in 1:col) {
        for (i in 1:(row - 1)) {
            for (j in (i + 1):row) {
                c = c + nvect[i, l] * nvect[j, l]
            }
        }
    }
    d = choose(N, 2) - a - b - c
    rad = (a + d)/choose(N, 2)
    ind <- c(rad, arj)
    names(ind) <- c("Rand Index (RI)", "Adjusted Rand Index (ARI)")
    ind
}
