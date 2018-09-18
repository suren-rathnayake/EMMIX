emmixfit1 <- function (dat, g, clust, distr, ncov, itmax, epsilon, initloop = 20)
{
    ndist <- switch(tolower(distr), mvn = 1, mvt = 2, msn = 3,
        mst = 4, 5)

    if (ndist > 4 | ncov < 1 | ncov > 5)
        stop("the model specified is not available yet")
    
    dat <- as.matrix(dat)
    
    n <- nrow(dat)
    p <- ncol(dat)
    
    if (n <= 20 * g)
        stop("sample size is too small")
    
    if (missing(clust) | is.null(clust))
        stop("initial clust must be given")
    
    clust <- unclass(as.ordered(clust))
    
    if (max(clust) != g)
        stop(paste("The levels of cluster should be g=", g))
    
    obj <- .C("emmixfit1", PACKAGE = "EMMIX", as.double(dat),
        as.integer(n), as.integer(p), as.integer(g), as.integer(ncov),
        as.integer(ndist), pro = double(g), mu = double(p * g),
        sigma = double(p * p * g), dof = double(g), delta = double(p *
            g), tau = double(n * g), double(n * g), double(n *
            g), double(n * g), double(n * g), sumtau = double(g),
        sumvt = double(g), sumzt = double(g), sumlnv = double(g),
        ewy = double(p * g), ewz = double(p * g), ewyy = double(p *
            p * g), loglik = double(1), lk = double(itmax), aic = double(1),
        bic = double(1), clust = as.integer(clust), error = integer(1),
        as.integer(itmax), as.double(epsilon), as.integer(initloop))[c(7:12,
        24:29)]
    
    lk <- obj$lk
    
    lk <- lk[lk != 0]
    
    list(distr = distr, error = obj$error, loglik = obj$loglik,
        bic = obj$bic, aic = obj$aic, pro = obj$pro, mu = array(obj$mu,
            c(p, g)), sigma = array(obj$sigma, c(p, p, g)), dof = obj$dof,
        delta = array(obj$delta, c(p, g)), clust = obj$clust,
        tau = array(obj$tau, c(n, g)), lk = lk)
}
