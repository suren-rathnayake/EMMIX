emmixfit2 <- function (dat, g, init, distr, ncov, itmax, epsilon)
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
    if (is.null(init) | missing(init))
        stop("init should be provided")
    pro <- init$pro
    mu <- init$mu
    sigma <- init$sigma
    dof <- init$dof
    delta <- init$delta
    obj <- .C("emmixfit2", PACKAGE = "EMMIX", as.double(dat),
        as.integer(n), as.integer(p), as.integer(g), as.integer(ncov),
        as.integer(ndist), pro = as.double(pro), mu = as.double(mu),
        sigma = as.double(sigma), dof = as.double(dof), delta = as.double(delta),
        tau = double(n * g), double(n * g), double(n * g), double(n *
            g), double(n * g), sumtau = double(g), sumvt = double(g),
        sumzt = double(g), sumlnv = double(g), loglik = double(1),
        lk = double(itmax), aic = double(1), bic = double(1),
        clust = integer(n), error = integer(1), as.integer(itmax),
        as.double(epsilon))[c(7:12, 21:26)]
    lk <- obj$lk
    lk <- lk[lk != 0]
    list(distr = distr, error = obj$error, loglik = obj$loglik,
        bic = obj$bic, aic = obj$aic, pro = obj$pro, mu = array(obj$mu,
            c(p, g)), sigma = array(obj$sigma, c(p, p, g)), dof = obj$dof,
        delta = array(obj$delta, c(p, g)), clust = obj$clust,
        tau = array(obj$tau, c(n, g)), lk = lk)
}
