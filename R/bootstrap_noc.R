bootstrap.noc <- function (x, g1, g2, distr, ncov, B = 99, replace = TRUE,
    itmax = 1000, epsilon = 1e-05, ...)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (g1 >= g2)
        stop("g1 should be less than g2")
    if (g1 < 1)
        stop("g1 should be greater than 0")
    counter <- 0
    vlk <- rep(0, g2 - g1 + 1)
    ret <- array(0, c(B, g2 - g1))
    dimnames(ret) <- list(1:B, paste(1 + g1:(g2 - 1), "vs", (g1:(g2 -
        1)), sep = " "))
    lk0 <- -Inf
    clust <- rep(1, n)
    for (g in g1:g2) {
        counter <- 0
        lk1 <- -Inf
        while (counter < 10) {
            if (g > 1)
                clust <- kmeans(x, g, nstart = 5)$cluster
            emobj <- emmixfit1(x, g, clust, distr, ncov, itmax,
                epsilon)
            if (emobj$error > 1)
                next
            if (emobj$loglik > lk1)
                lk1 <- emobj$loglik
            counter = counter + 1
        }
        # dput(emobj, paste("ReturnOf_g_", g, ".ret", sep = ""))

        counter <- 0
        lk0 <- lk1
        vlk[g - g1 + 1] <- lk0
        if (g < g2) {
            for (i in 1:(2 * B)) {
                if (replace)
                  dat <- x[sample(1:n, n, replace = TRUE), ]
                else dat <- rdemmix2(n, distr, emobj$pro,
                  emobj$mu, emobj$sigma, emobj$dof, emobj$delta)
                if (is.null(dat))
                  stop("I can not generate the data!")
                obj <- emmixfit2(dat, g, emobj, distr, ncov,
                  itmax, epsilon)
                if (obj$error > 1)
                  next
                ii <- 0
                lk2 <- -Inf
                while (ii < 10) {
                  clust <- kmeans(dat, g + 1, nstart = 5)$cluster
                  obj2 <- emmixfit1(dat, g + 1, clust, distr,
                    ncov, itmax, epsilon)
                  ii <- ii + 1
                  if (obj2$error > 1)
                    next
                  if (obj2$loglik > lk2)
                    lk2 <- obj2$loglik
                }
                counter <- counter + 1
                ret[counter, g - g1 + 1] <- -2 * (obj$loglik -
                  lk2)
                if (counter >= B)
                  break
            }
        }
    }
    pvalue <- rep(0, g2 - g1)
    for (i in 1:(g2 - g1)) {
        pvalue[i] <- sum(ret[, i] < 2 * (vlk[i + 1] - vlk[i]))/B
    }
    list(lrts = ret[], log_lk = vlk, pvalue = pvalue)
}