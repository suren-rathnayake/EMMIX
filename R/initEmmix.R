initEmmix <- function (dat, g, clust, distr, ncov, maxloop = 20)
{
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3,
        mst = 4, 5)
    if (ndist > 4 | ncov < 1 | ncov > 5)
        stop("the model specified is not available yet")
    dat <- as.matrix(dat)
    n <- nrow(dat)
    p <- ncol(dat)
    clust <- unclass(as.ordered(clust))
    obj <- .C("initfit_", PACKAGE = "EMMIX", as.double(dat),
        as.integer(n), as.integer(p), as.integer(g), as.integer(ncov),
        as.integer(ndist), pro = double(g), mu = double(p * g),
        sigma = double(p * p * g), dof = double(g), delta = double(p *
            g), tau = double(n * g), double(n * g), double(n *
            g), double(n * g), double(n * g), sumtau = double(g),
        sumvt = double(g), sumzt = double(g), sumlnv = double(g),
        ewy = double(p * g), ewz = double(p * g), ewyy = double(p *
            p * g), loglik = double(1), as.integer(clust), error = integer(1),
        as.integer(maxloop))[c(7:11, 24, 26)]
    error <- obj$error
    ret <- NULL
    if (error == 0) {
        ret <- list(distr = distr, error = error, loglik = obj$loglik,
            pro = obj$pro, mu = array(obj$mu, c(p, g)), sigma = array(obj$sigma,
                c(p, p, g)), dof = obj$dof, delta = array(obj$delta,
                c(p, g)))
    }
    else warning("error:", error)
    ret
}
