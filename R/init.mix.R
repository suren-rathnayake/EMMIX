init.mix <- function (dat, g, distr, ncov, nkmeans, nrandom, nhclust, maxloop = 20)
{
    found <- list()
    found$loglik <- -Inf
    n <- nrow(dat)
    clust <- rep(1, n)
    mclust <- NULL
    if (g > 1) {
        if (nkmeans > 0) {
            for (i in 1:nkmeans) {
                clust <- kmeans(dat, g, nstart = 5)$cluster
                if (min(table(clust)) < 10)
                  next
                obj <- try(initEmmix(dat, g, clust, distr, ncov,
                  maxloop))
                if (length(obj) != 8 | obj$error)
                  next
                if (obj$loglik > found$loglik) {
                  found <- obj
                  mclust <- clust
                }
            }
            if (is.null(mclust))
                nrandom <- 10
        }
        if (nrandom > 0)
            for (i in 1:nrandom) {
                clust <- sample(1:g, n, replace = TRUE)
                obj <- try(initEmmix(dat, g, clust, distr, ncov,
                  maxloop))
                if (length(obj) != 8 | obj$error != 0)
                  next
                if (obj$loglik > found$loglik) {
                  found <- obj
                  mclust <- clust
                }
            }
        methods <- c("complete")
        if (nhclust) {
            dd <- as.dist((1 - cor(t(dat)))/2)
            for (j in methods) {
                clust <- cutree(hclust(dd, j), k = g)
                obj <- try(initEmmix(dat, g, clust, distr, ncov,
                  maxloop))
                if (length(obj) != 8 | obj$error != 0)
                  next
                if (obj$loglik > found$loglik) {
                  found <- obj
                  mclust <- clust
                }
            }
        }
    }
    else {
        obj <- try(initEmmix(dat, g, clust, distr, ncov, maxloop))
        if (length(obj) == 8) {
            found <- obj
            mclust <- clust
        }
    }
    if (is.null(mclust)) {
        found <- NULL
        warning("failed to find a initial values!")
    }
    found
}
