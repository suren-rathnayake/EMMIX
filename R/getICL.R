getICL <- function (x, n, p, g, distr, ncov, pro, mu, sigma, dof, delta,
    clust)
{
    x <- as.matrix(x)
    loglik <- 0
    nc = 0
    nn <- sum(clust > 0)
    lnden <- (as.matrix((ddmix(x, n, p, g, distr, mu, sigma,
        dof, delta))))
    for (h in 1:g) loglik = loglik + sum(ifelse(clust == h, log(pro[h]) +
        lnden[, h], 0))
    nc <- switch(tolower(ncov), `1` = (g - 1) + g * p + p * (1 +
        p)/2, `2` = (g - 1) + g * p + p, `3` = (g - 1) + g *
        (p + p * (1 + p)/2), `4` = (g - 1) + g * (p + p), `5` = (g -
        1) + g * (p + 1))
    nc <- switch(tolower(distr), mvn = nc, mvt = nc + g, msn = nc +
        g * p, mst = nc + g * p + g)
    ICL = loglik - nc/2 * log(nn)
    list(ICL = ICL)
}
