rdmvn <- function (n, mean, cov)
{
    cov <- as.matrix(cov)
    p <- nrow(cov)

    if (nrow(cov) != ncol(cov)) {
        stop("cov must be a square matrix")
    }
    if (length(mean) != nrow(cov)) {
        stop("mean and cov have non-conforming size")
    }
    rmvnorm(n, mean = mean, sigma = cov, method = "chol")
}


rdmvt <- function (n, mean, cov, nu)
{
    cov <- as.matrix(cov)
    p <- nrow(cov)
    u <- rgamma(n, nu/2, nu/2)
    t(t(rdmvn(n, mean = rep(0, p), cov = cov)/sqrt(u)) + mean)
}


rdemmix <- function (nvect, p, g, distr, mu, sigma, dof = NULL, delta = NULL)
{
    if (length(c(nvect)) != g)
        stop("nvect should be a vector")
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3,
        mst = 4, 5)
    if (ndist > 4)
        stop("the model specified is not available yet")
    if (is.null(dof))
        dof <- rep(4, g)
    if (is.null(delta))
        delta <- array(0, c(p, g))
    if (length(c(mu)) != (p * g))
        stop(paste("mu should be a ", p, "by", g, "matrix!"))
    if (length(c(sigma)) != (p * p * g))
        stop(paste("sigma should be a ", p, "by", p, "by", g,
            " array!"))
    if (length(c(dof)) != g)
        stop(paste("dof should be a ", g, " vector!"))
    if (length(c(delta)) != (p * g))
        stop(paste("delta should be a ", p, "by", g, " array!"))
    mu = array(mu, c(p, g))
    sigma = array(sigma, c(p, p, g))
    delta = array(delta, c(p, g))
    dat <- array(0, c(10, p))
    mvrand <- function(n, p, ndist, mean, cov, nu, del) {
        switch(ndist, `1` = rdmvn(n, mean = mean, cov = cov),
            `2` = rdmvt(n, p, mean = mean, cov = cov, nu = nu),
            `3` = rdmsn(n, p, mean = mean, cov = cov, del = del),
            `4` = rdmst(n, p, mean = mean, cov = cov, nu = nu,
                del = del))
    }
    if (g >= 1)
        for (h in 1:g) {
            if (nvect[h] > 0)
                dat <- rbind(dat, mvrand(nvect[h], p, ndist,
                  mu[, h], sigma[, , h], dof[h], delta[, h]))
        }
    dat[-(1:10), ]
}


rdemmix2 <- function (n, p, g, distr, pro, mu, sigma, dof = NULL, delta = NULL)
{
    n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
    nn <- n0
    if (length(nn) < g) {
        nn <- rep(0, g)
        for (i in as.numeric(names(n0))) nn[i] <- n0[paste(i)]
    }
    names(nn) <- NULL
    rdemmix(nn, p, g, distr, mu, sigma, dof, delta)
}



rdemmix3 <- function (n, p, g, distr, pro, mu, sigma, dof = NULL, delta = NULL)
{
    if (length(pro) != g)
        stop(paste("pro should be a ", g, " vector!"))
    n0 <- table(sample(1:g, n, replace = TRUE, prob = pro))
    nn <- n0
    if (length(nn) < g) {
        nn <- rep(0, g)
        for (i in as.numeric(names(n0))) nn[i] <- n0[paste(i)]
    }
    names(nn) <- NULL
    dat <- rdemmix(nn, p, g, distr, mu, sigma, dof, delta)
    list(data = dat, cluster = rep(1:g, nn))
}

ddmvn <- function (dat, mean, cov)
{
    p <- nrow(cov)
    exp(ddmix(dat, 1, "mvn", mean, cov, 0, rep(0, p)))
}

ddmvt <- function (dat, mean, cov, nu)
{
  p <- nrow(cov)
  exp(ddmix(dat, 1, "mvt", mean, cov, nu, rep(0, p)))
}

# ddmix <-
# function (dat, g, distr, mu, sigma, dof = NULL, delta = NULL)
# {
#     n <- nrow(x)
#     p <- ncol(x)
#
#     if (is.null(dof))
#         dof <- rep(4, g)
#     if (is.null(delta))
#         delta <- array(0, c(p, g))
#     ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3,
#         mst = 4, 5)
#     if (ndist > 4)
#         stop("the model specified is not available yet")
#     dat <- as.matrix(dat)
#     if (n == 1 & (ncol(dat) == 1))
#         dat <- t(dat)
#     if (nrow(dat) != n | ncol(dat) != p)
#         stop("dat does not match n and p.")
#     if (length(c(mu)) != (p * g))
#         stop(paste("mu should be a ", p, "by", g, "matrix!"))
#     if (length(c(sigma)) != (p * p * g))
#         stop(paste("sigma should be a ", p, "by", p, "by", g,
#             " array!"))
#     if (length(c(dof)) != g)
#         stop(paste("dof should be a ", g, " vector!"))
#     if (length(c(delta)) != (p * g))
#         stop(paste("delta should be a ", p, "by", g, " array!"))
#     obj <- .C("ddmix", PACKAGE = "EMMIX", as.double(dat), as.integer(n),
#         as.integer(p), as.integer(g), as.integer(ndist), as.double(mu),
#         as.double(sigma), as.double(dof), as.double(delta), den = double(n *
#             g), error = integer(1))[10:11]
#     if (obj$error)
#         stop("error")
#     (matrix(obj$den, ncol = g))
# }

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
